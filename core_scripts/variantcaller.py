import re
import os
import copy
import pandas as pd
from tools.set_env import Environment
from tools.logger import Logger
from tools.find_refstartpos import get_ref_startpos
from tools.names_generator import create_vcf_name_from_pair_and_ref
from constants import FASTA_CHRS_PATH, REF_NAMES, PRIMARY_PATH, \
    ANCESTRAL_COLNAMES, LOG_PATH
from pyliftover import LiftOver


# wrapper to extract SNP's
def variant_call_pairs(pairs_list, chr_to_select, run_id=None):
    caller_log = Logger(reset=True, source_name='variant_call_pairs',
                        run_id=run_id)
    chr_names = chr_to_select
    processed_files = []
    for pair in pairs_list:
        caller_log.tick(timer_name='pair_timer')
        caller_log.message(msg='Started pair ' + str(pair.string()))
        v_df1, v_df2 = [pd.DataFrame(), pd.DataFrame()]
        for chrname in chr_names:
            caller_log.tick(timer_name='chr_timer')
            caller_log.message(msg='Calling variants for ' + chrname)

            if pair.reference not in REF_NAMES:
                chrv1, chrv2 = chr_variant_call(chrname,
                                                sp_name_1=pair.species_1,
                                                sp_name_2=pair.species_2,
                                                ref_name=pair.reference,
                                                run_id=run_id)
            else:
                chrv1, chrv2 = chr_variant_call_reconstructed(chrname,
                                                              sp_name_1=pair.species_1,
                                                              sp_name_2=pair.species_2,
                                                              ref_name=pair.reference,
                                                              run_id=run_id)

            caller_log.message(msg="finished chr " + chrname, print_time=True,
                               timer_name='chr_timer')
            v_df1 = pd.concat([chrv1, v_df1], ignore_index=True)
            v_df2 = pd.concat([chrv2, v_df2], ignore_index=True)
        caller_log.message(msg='Ended pair ' + str(pair.string()),
                           print_time=True,
                           timer_name='pair_timer')
        sp1_filename = df_to_vcf(v_df1, sp_name=pair.species_1,
                                 sp_pair=pair.species_2,
                                 refname=pair.reference,
                                 run_id=run_id)
        sp2_filename = df_to_vcf(v_df2, sp_name=pair.species_2,
                                 sp_pair=pair.species_1,
                                 refname=pair.reference,
                                 run_id=run_id)
        processed_files.append(sp1_filename)
        processed_files.append(sp2_filename)
    caller_log.print_end()
    return processed_files


def chr_variant_call(chr_name, sp_name_1="", sp_name_2="", ref_name='',
                     run_id=None):
    outputpath = os.path.join(FASTA_CHRS_PATH, chr_name)
    variant_call_logger = Logger(source_name='variant_call', run_id=run_id)
    #                             msg='variants counting for ' + chr_name + ' ' + sp_name_1 + ' ' + sp_name_2)
    sp_names = [ref_name, sp_name_1, sp_name_2]
    f = [open(os.path.join(outputpath, name + ".fasta"), 'r') for name in
         sp_names]
    files = dict(zip(sp_names, f))
    docpos = 0
    alt_name = ''
    context = ''
    [ALT, REF, REFPOS, INFO, CHR, CONTEXT] = [[], [], [], [], [], []]
    for name in files.keys():
        files[name].readline()
    [sp_1_prev, sp_2_prev, ref_prev] = [files[sp_name_1].read(1).upper(),
                                        files[sp_name_2].read(1).upper(),
                                        files[ref_name].read(1).upper()]
    [sp_1_nuc, sp_2_nuc, ref_nuc] = [files[sp_name_1].read(1).upper(),
                                     files[sp_name_2].read(1).upper(),
                                     files[ref_name].read(1).upper()]
    refpos = get_ref_startpos(chr_name)
    docpos += 1
    refpos += 1

    while True:
        [sp_1_next, sp_2_next, ref_next] = [files[sp_name_1].read(1).upper(),
                                            files[sp_name_2].read(1).upper(),
                                            files[ref_name].read(1).upper()]
        refpos += 1
        alt = 'Q'
        docpos += 1
        if bool(re.match(r'[ATGCatgc]', sp_1_nuc)) and \
                bool(re.match(r'[ATGCatgc]', sp_2_nuc)) and \
                bool(re.match(r'[ATGCatgc]', ref_nuc)):
            if sp_1_nuc != sp_2_nuc:
                if not (sp_1_nuc != ref_nuc and sp_2_nuc != ref_nuc):
                    context_ref = ref_prev + '.' + ref_next
                    if sp_1_nuc == ref_nuc:
                        alt = sp_2_nuc
                        alt_name = sp_name_2
                        context = sp_2_prev + '.' + sp_2_next
                        context_alt = sp_1_prev + '.' + sp_1_next
                    elif sp_2_nuc == ref_nuc:
                        alt = sp_1_nuc
                        alt_name = sp_name_1
                        context = sp_1_prev + '.' + sp_1_next
                        context_alt = sp_2_prev + '.' + sp_2_next
                    else:
                        variant_call_logger.message(
                            msg='Error when comparing ')
                    if bool(re.match(r'[ATGCatgc].[ATGCatgc]', context)):
                        info = str(docpos)
                        CONTEXT.append(context.upper())
                        CHR.append(chr_name + '_' + alt_name)
                        ALT.append(alt)
                        REF.append(ref_nuc)
                        REFPOS.append(refpos)
                        if not (bool(re.match(r'[ATGCatgc].[ATGCatgc]', \
                                              context_alt))) or \
                                not (bool(re.match(r'[ATGCatgc].[ATGCatgc]', \
                                                   context_ref))):
                            info += '_badcontext'
                        INFO.append(info)
        [sp_1_prev, sp_2_prev, ref_prev] = [sp_1_nuc, sp_2_nuc, ref_nuc]
        [sp_1_nuc, sp_2_nuc, ref_nuc] = [sp_1_next, sp_2_next, ref_next]

        if not sp_1_nuc or not sp_2_nuc:
            break
    SUBS = [nref.upper() + '>' + nmut.upper() for nref, nmut in zip(REF, ALT)]
    for index, sbs in enumerate(SUBS):
        if sbs == 'A>T':
            SUBS[index] = 'T>A'
            INFO[index] += (' - A>T')
        if sbs == 'G>C':
            SUBS[index] = 'C>G'
            INFO[index] += (' - G>C')
        if sbs == 'G>A':
            SUBS[index] = 'C>T'
            INFO[index] += (' - G>A')
        if sbs == 'A>G':
            SUBS[index] = 'T>C'
            INFO[index] += (' - A>G')
        if sbs == 'A>C':
            SUBS[index] = 'T>G'
            INFO[index] += (' - A>C')
        if sbs == 'G>T':
            SUBS[index] = 'C>A'
            INFO[index] += (' - G>T')

    dictdata = {'Chr': CHR, 'SUBS': SUBS, 'Alt': ALT, 'Ref': REF,
                'Refpos': REFPOS, 'Context': CONTEXT, 'Info': INFO}
    chrdata = pd.DataFrame(dictdata, columns=dictdata.keys())
    data_sp_1 = pd.DataFrame(
        chrdata.loc[chrdata['Chr'] == (chr_name + '_' + sp_name_1)])
    data_sp_2 = pd.DataFrame(
        chrdata.loc[chrdata['Chr'] == (chr_name + '_' + sp_name_2)])

    return data_sp_1, data_sp_2


def chr_variant_call_reconstructed(chr_name, sp_name_1="", sp_name_2="",
                                   ref_name='',
                                   run_id=None):
    path_to_liftovers = os.path.join(PRIMARY_PATH, 'data', 'liftovers')
    outputpath = os.path.join(FASTA_CHRS_PATH, chr_name)
    ancestral_alleles_path = os.path.join(PRIMARY_PATH, 'data',
                                          'ancestral_vcfs')
    variant_call_logger = Logger(source_name='variant_call anc', run_id=run_id)
    alt_name, context = ['', '']
    [ALT, REF, REFPOS, INFO, CHR, CONTEXT] = [[], [], [], [], [], []]
    [not_liftovered_bases, bases_not_matched] = [[], []]
    [not_liftovered_pos_counter, base_counter, bases_matched_counter,
     hg18_coord] \
        = [0, 0, 0, 0]

    # loading species fasta
    sp_names = [sp_name_1, sp_name_2]
    f = [open(os.path.join(outputpath, name + ".fasta"), 'r') for name in
         sp_names]
    files = dict(zip(sp_names, f))
    for name in files.keys():
        files[name].readline()

    # configuring liftover
    lo18to38 = LiftOver(
        os.path.join(path_to_liftovers, 'hg18ToHg38.over.chain.gz'))

    # loading ancestral allels
    ref_column_index = \
        [ind for ind, colname in enumerate(ANCESTRAL_COLNAMES) if
         colname == ref_name][0]

    ref_df = pd.read_csv(os.path.join(
        ancestral_alleles_path, chr_name + '.calls.txt'),
        skiprows=1, usecols=[1, ref_column_index], sep='\t')
    ref_df.columns = ['position', ref_name]
    # ref_df = ref_df.set_index(['position'],drop=False)

    ref_df['hg38_pos'] = pd.Series(
        [lo18to38.convert_coordinate(chr_name, x)[0][1]
         if lo18to38.convert_coordinate(chr_name, x) != []
         else -1
         for x in ref_df['position']]).values
    hg38_coord_min = get_ref_startpos(chr_name)
    # sorting by liftovered coords and removind indefinite states
    ref_df = ref_df.sort_values(by=['hg38_pos'], ascending=True)
    ref_df = ref_df.drop(ref_df[ref_df['hg38_pos'] == -1].index)
    ref_df = ref_df.loc[
        ref_df[ref_name].str.contains('^[ATGCatgc]$', regex=True)]

    ref_df = ref_df.drop_duplicates(['hg38_pos'])
    ref_pos_min = ref_df['hg38_pos'].min()
    ref_df = ref_df.drop(ref_df[ref_df['hg38_pos'] <= hg38_coord_min].index)

    ref_df = ref_df.reset_index(drop=True)
    print('ref_df_shape ', ref_df.shape)

    ref_df_indexer = 0
    max_ref_pos = ref_df['hg38_pos'].values[-1]
    hg38_coord = get_ref_startpos(chr_name)



    # start looping in file
    [sp_1_prev, sp_2_prev] = [files[sp_name_1].read(1).upper(),
                              files[sp_name_2].read(1).upper()]
    [sp_1_nuc, sp_2_nuc] = [files[sp_name_1].read(1).upper(),
                            files[sp_name_2].read(1).upper()]
    ref_nuc = ref_df.loc[ref_df_indexer, ref_name]
    ref_pos = ref_df.loc[ref_df_indexer, 'hg38_pos']

    hg38_coord += 1
    while True:
        [sp_1_next, sp_2_next] = [files[sp_name_1].read(1).upper(),
                                  files[sp_name_2].read(1).upper()]
        hg38_coord += 1
        base_counter += 1
        alt = 'undef_n'
        # comapre only for coords that present in mrca
        if hg38_coord > ref_pos:
            variant_call_logger.message(
                msg='Error when comparing ref_df index {df_index} ref_pos {ref_pos} hg38coord {hg38_coord} '
                    .format(df_index=ref_df_indexer,
                            hg38_coord=hg38_coord,
                            ref_pos=ref_pos))
        if hg38_coord == ref_pos:
            if bool(re.match(r'[ATGCatgc]', sp_1_nuc)) and \
                    bool(re.match(r'[ATGCatgc]', sp_2_nuc)):
                if sp_1_nuc != sp_2_nuc:
                    bases_matched_counter += 1
                    if not (sp_1_nuc != ref_nuc and sp_2_nuc != ref_nuc):
                        if sp_1_nuc == ref_nuc:
                            alt = sp_2_nuc
                            alt_name = sp_name_2
                            context = sp_2_prev + '.' + sp_2_next
                            context_alt = sp_1_prev + '.' + sp_1_next
                        elif sp_2_nuc == ref_nuc:
                            alt = sp_1_nuc
                            alt_name = sp_name_1
                            context = sp_1_prev + '.' + sp_1_next
                            context_alt = sp_2_prev + '.' + sp_2_next
                        else:
                            variant_call_logger.message(
                                msg='Error when comparing {chr_name}:{hg38_coord} position'
                                    .format(chr_name=chr_name,
                                            hg38_coord=hg38_coord))
                        if bool(re.match(r'[ATGCatgc].[ATGCatgc]', context)):
                            info = str(base_counter)
                            CONTEXT.append(context.upper())
                            CHR.append(chr_name + '_' + alt_name)
                            ALT.append(alt)
                            REF.append(ref_nuc)
                            REFPOS.append(hg38_coord)
                            if not (bool(re.match(r'[ATGCatgc].[ATGCatgc]', \
                                                  context_alt))):
                                info += '_badcontext'
                            INFO.append(info)
            ref_df_indexer += 1
            ref_nuc = ref_df.loc[ref_df_indexer, ref_name]
            ref_pos = ref_df.loc[ref_df_indexer, 'hg38_pos']
        # next step indide the file
        [sp_1_prev, sp_2_prev] = [sp_1_nuc, sp_2_nuc]
        [sp_1_nuc, sp_2_nuc] = [sp_1_next, sp_2_next]

        if hg38_coord % 20000000 == 0:
            variant_call_logger.message(
                'passed {chr_name}: {hg38_coord} ' 
                ' {matched_bases} bases matched'.format(
                    chr_name=chr_name, hg38_coord=hg38_coord,
                matched_bases=bases_matched_counter))

        if not sp_1_nuc or not sp_2_nuc or ref_pos == max_ref_pos:
            print('stopped at ', hg38_coord)
            break

    SUBS = [nref.upper() + '>' + nmut.upper() for nref, nmut in zip(REF, ALT)]
    dictdata = {'Chr': CHR, 'SUBS': SUBS, 'Alt': ALT, 'Ref': REF,
                'Refpos': REFPOS, 'Context': CONTEXT, 'Info': INFO}
    chrdata = pd.DataFrame(dictdata, columns=dictdata.keys())

    if len(SUBS) < 6:
        variant_call_logger.message(
            'No varinats found for {chr_name}'.format(chr_name=chr_name))
        return chrdata, chrdata

    data_sp_1 = pd.DataFrame(
        chrdata.loc[chrdata['Chr'] == (chr_name + '_' + sp_name_1)])
    data_sp_2 = pd.DataFrame(
        chrdata.loc[chrdata['Chr'] == (chr_name + '_' + sp_name_2)])

    variant_call_logger.message('there are {matched} matched bases'
                                .format(matched=bases_matched_counter))

    return data_sp_1, data_sp_2


def df_to_vcf(result, sp_name, sp_pair, refname, run_id=None):
    final_pair_path = os.path.join(Environment().variant_path, 'raw_variants')
    if not os.path.exists(final_pair_path):
        os.makedirs(final_pair_path)
    filename = create_vcf_name_from_pair_and_ref(sp_name=sp_name,
                                                 sp_pair=sp_pair,
                                                 refname=refname)
    header = '##fileformat=VCFv4.1 \n'
    header += (
            '##contig=<ID=' + filename + ', length=0,assembly=hg38> \n##reference=' + refname + '\n')
    header += '#CHROM POS ID REF ALT QUAL FILTER INFO\n'
    f = open(filename, 'w')
    f.write(header)
    f.close()

    # del result['Info']
    result["CHROM"] = result['Chr'].replace(r'_\w+', "", inplace=False,
                                            regex=True)
    result['Chr'] = result["Chr"] + "_" + result["SUBS"] + "_" + result[
        "Context"]
    result['QUAL'] = '.'
    result['FILTER'] = '.'
    # result['INFO'] = '.'
    del result['Context']
    result.rename(index=str,
                  columns={"Chr": "ID", "Refpos": "POS", "Ref": "REF",
                           "Alt": "ALT", "Info": 'INFO'}, inplace=True)
    result = result[
        ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"]]
    f = open(filename, 'r')
    HH = []
    for line in f:
        HH.append(line.strip())
    # print(HH)
    f.close()
    with open(filename, 'a') as f:
        result.to_csv(f, header=False, sep='\t', index=False, mode='a')
    logher = Logger(source_name='to_vcf_convertion',
                    msg='finished vcf for {}'.format(sp_name), run_id=run_id)
    return filename
