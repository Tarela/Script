import os
import yaml
from operator import itemgetter, eq
import Levenshtein
import logging
# create logger
logger = logging.getLogger("get_barcodes")
logger.setLevel(logging.DEBUG)
# create console handler and set level to debug
fh = logging.FileHandler('2interesting_log_get_barcodes.log') # lucky number 7!!!!! Give me power!!!
fh.setLevel(logging.DEBUG)
# create formatter
formatter = logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s")
# add formatter to ch
fh.setFormatter(formatter)
# add ch to logger
logger.addHandler(fh)

# for mapping add @
def add_tag_to_fastq(fastq_file):
    modified_fastq = 'modified_' + os.path.basename(fastq_file)
    modified_fastq_path = os.path.join(os.path.dirname(fastq_file), modified_fastq)
    cnt = 0
    with open(fastq_file, 'r') as fopen:
        with open(modified_fastq_path, 'w') as mopen:
            for line in fopen:
                cnt += 1
                if cnt % 4 == 1:
                    mopen.write('@' + line)
                else:
                    mopen.write(line)

# get barcodes from indextable.txt
def get_barcodes_templates(index_table):
    barcodes_index_list = []
    with open(index_table, 'r') as fopen:
        line = fopen.readline() # skip the first line
        line = fopen.readline()
        while line:
            barcodes_index_list.append(line.strip().split('\t')[0])
            line = fopen.readline()
    return barcodes_index_list

# get barcodes from fastq
def get_observed_barcodes(fastq_file):
    count = 0
    barcodes_dict = {}
    barcodes_list = []
    with open(fastq_file, 'r') as fopen:
        line = fopen.readline().strip()
        while line:
            # logger.info(line)
            count += 1
            barcode = line.split(':')[0] # get barcode
            barcodes_list.append(barcode)
            if barcode in barcodes_dict:
                barcodes_dict[barcode] += 1
            else:
                barcodes_dict[barcode] = 1
            # read three lines
            fopen.readline() # sequence
            fopen.readline() # tag
            fopen.readline() # phred score
            line = fopen.readline().strip()
    return barcodes_dict, barcodes_list


# get barcodes from readcount
def get_readcount_barcodes(readcounts_file):
    rc_barcodes_list = []
    with open(readcounts_file, 'r') as fopen:
        line = fopen.readline() # skip the first line
        line = fopen.readline()
        while line:
            readcount_line = line.strip().split('\t')
            rc_barcodes_list.append((readcount_line[0], readcount_line[2])) # (seq, totalreads)
            line = fopen.readline()
    return rc_barcodes_list

# split barcodes from files
def split_barcodes(indextable_file):
    barcodes_dict = {'trans1': [], 'pcr1': [], 'pcr2':[], 'trans2': []}
    with open(indextable_file, 'r') as fopen:
        line = fopen.readline() # skip the first line
        line = fopen.readline()
        while line:
            barcodes = line.strip().split('\t')[0]
            if barcodes[0:8] not in barcodes_dict['trans1']:
                barcodes_dict['trans1'].append(barcodes[0:8])
            if barcodes[8:16] not in barcodes_dict['pcr1']:
                barcodes_dict['pcr1'].append(barcodes[8:16])
            if barcodes[16:24] not in barcodes_dict['pcr2']:
                barcodes_dict['pcr2'].append(barcodes[16:24])
            if barcodes[24:32] not in barcodes_dict['trans2']:
                barcodes_dict['trans2'].append(barcodes[24:32])
            line = fopen.readline()

    return barcodes_dict

# add reverse into account
def amplify_barcodes_dict(bc_dict):
    reverse_dict = {'trans1': 'trans2', 'trans2': 'trans1', 'pcr1': 'pcr2', 'pcr2': 'pcr1'}

    reverse_bc_dict = {}
    for bc_type, bc_list in bc_dict.items():
        reverse_list = []
        for bc in bc_list:
            bc = bc[::-1]
            reverse_list.append(bc.replace('A', 't').replace('T','a').replace('G', 'c').replace('C', 'g').upper())
        reverse_bc_dict[reverse_dict[bc_type]] = reverse_list
    
    for bc_type, bc_list in bc_dict.items():
        bc_list.extend(reverse_bc_dict[bc_type])
    

# according to the shendure lab rule
def generate_indexconversion_table(barcodes_dict, barcodes_fastq_list, ic_list):
    # amplify_barcodes_dict(barcodes_dict)
    logger.info(barcodes_dict)
    logger.info('trans1: %s' % len(barcodes_dict['trans1']))
    logger.info('trans2: %s' % len(barcodes_dict['trans2']))
    logger.info('pcr1: %s' % len(barcodes_dict['pcr1']))
    logger.info('pcr2: %s' % len(barcodes_dict['pcr2']))

    trans1 = barcodes_dict['trans1']
    trans2 = barcodes_dict['trans2']
    pcr1 = barcodes_dict['pcr1']
    pcr2 = barcodes_dict['pcr2']


    # shit everywhere!
    # condition_dict = {'no_rule': [], 'less_than_2': [], 'less_than_2_2': [], 'less_than_3': [], 'less_than_3_and_1': [], 'less_than_3_and_2': [], 'paper': []}
    # condition_dict = {'paper': []}
    condition_dict = {'less_than_3_and_2': [], 'paper': []}
    for condition in condition_dict.keys(): # shit
        combo = []
        for observed_barcode in barcodes_fastq_list:
            ob_1 = observed_barcode[0:8]
            ob_2 = observed_barcode[8:16]
            ob_3 = observed_barcode[16:24]
            ob_4 = observed_barcode[24:32]

            # if reverse complement
            # ob_3 = observed_barcode[24:32][::-1].replace('A', 't').replace('T','a').replace('G', 'c').replace('C', 'g').upper()
            # ob_4 = observed_barcode[16:24][::-1].replace('A', 't').replace('T','a').replace('G', 'c').replace('C', 'g').upper()

        
            # here comes some shit!
            new_barcode = map_ob_with_templates(ob_1, trans1, condition)
            new_barcode += map_ob_with_templates(ob_2, pcr1, condition)
            new_barcode += map_ob_with_templates(ob_3, pcr2, condition)
            new_barcode += map_ob_with_templates(ob_4, trans2, condition)
        
            # map to indexconversion table
            combo.append((observed_barcode, new_barcode))
        cmp_gen_with_file(combo, ic_list, condition) # fucking return [map, nomap]
    return combo

# shit comes
def map_ob_with_templates(ob, template_list, condition): 
    if ob in template_list:
        return ob

    # fly-atac
    if condition == 'paper':
        winner = '_CTF' + '_'*(len(ob)-4)
        winner_ed = 8
        runnerup_ed = 8
        for tmp in template_list:
            curred = Levenshtein.distance(ob,tmp)
            if curred <= winner_ed:
                runnerup_ed = winner_ed
                winner = tmp
                winner_ed = curred
        if winner_ed > 3:
             winner = '_CTF' + '_'*(len(ob)-4)
        if runnerup_ed - winner_ed < 2:
            winner = '_AMBIG' + '_'*(len(ob)-6)
        return winner

    if condition == 'no_rule':
        return '_CTF____'
    
    dist_list = []
    for tmp in template_list:
        dist_list.append((tmp,Levenshtein.distance(ob, tmp)))

    dist_list = sorted(dist_list, key=itemgetter(1))
    # logger.info(ob)
    # logger.info(dist_list)
    # whether it fits the match 
    # shity shit
    if condition == 'less_than_3': 
        if dist_list[0][1] > 3: # which to choose as first one
            return '_CTF____'
        else:
            return dist_list[0][0]
    
    if condition == 'less_than_3_and_1': 
        if dist_list[0][1] > 3:
            return '_CTF____'
        else:
            if (dist_list[1][1] - dist_list[0][1]) < 1:
                return '_AMBIG__'
            if (dist_list[1][1] - dist_list[0][1]) >= 1:
                return dist_list[0][0]

    if condition == 'less_than_3_and_2': # what paper asks
        cur_ob = dist_list[0][0]
        if dist_list[0][1] > 3:
            cur_ob = '_CTF____'
        # interesting part
        if (dist_list[1][1] - dist_list[0][1]) < 2:
            cur_ob = '_AMBIG__'
        return cur_ob
        # else:
            # if (dist_list[1][1] - dist_list[0][1]) < 2:
            #    return '_AMBIG__'
            # if (dist_list[1][1] - dist_list[0][1]) >= 2:
            #     return dist_list[0][0]

    if condition == 'less_than_2': 
        if dist_list[0][1] > 2:
            return '_CTF____'
        else:
            if (dist_list[1][1] - dist_list[0][1]) < 2:
                return '_AMBIG__'
            if (dist_list[1][1] - dist_list[0][1]) >= 2:
                return dist_list[0][0]

    if condition == 'less_than_2_2':
        if dist_list[0][1] > 2:
            return '_CTF____'
        else:
            if (dist_list[1][1] - dist_list[0][1]) < 1:
                return '_AMBIG__'
            if (dist_list[1][1] - dist_list[0][1]) >= 1:
                return dist_list[0][0]



def cmp_gen_with_file(combo, ic_list, condition):
    # what if less than 3 and additional < 0, choose whatever

    logger.info(condition)
    logger.info('combo: %s' %  len(combo))
    logger.info('ic_list: %s' %  len(ic_list))
    logger.info('compare result: ')

    # tmp_combo = combo
    # tmp_ic = ic_list

    # for cmb in tmp_combo:
        # if cmb in ic_list:
            # logger.info(cmb)
            # tmp_combo.remove(cmb)
            # ic_list.remove(cmb)  # shit to shit
    ic_set = set(ic_list)
    logger.info('ic_set: %s ' % len(ic_set))
    not_map = [x for x in combo if x not in ic_set]
    logger.info(len(not_map))
    logger.info(not_map)
    return [len(not_map), len(combo)-len(not_map)]


def get_index_conversion(ic_file):
    ic_list = []
    with open(ic_file, 'r') as fopen:
        line = fopen.readline()
        line = fopen.readline()
        while line:
            number, origin, after = line.strip().split('\t')
            ic_list.append((origin, after))
            line = fopen.readline()
    return ic_list

def get_ed_matrix(barcodes_dict):
    for bc_type, bc_list in barcodes_dict.items():
        logger.info('current type: %s' % bc_type)
        matrix = [[0 for x in range(len(bc_list))] for y in range(len(bc_list))]
        for i, bc_x in enumerate(bc_list):
            log_str = ''
            for j, bc_y in enumerate(bc_list):
                matrix[i][j] = Levenshtein.distance(bc_x, bc_y)
                log_str += str(matrix[i][j]) + ' '
            logger.info(log_str)


if __name__ == '__main__':
    # HL_path = '/nv/vol190/zanglab/wm9tr/projects/scRNA_seq_beginner/barcode/data/GM12878vsHL'
    # fastq_file = 'SRR1947693_2.fastq'
    # matrix_file = 'GSM1647123_GM12878vsHL.dhsmatrix.txt'
    # indextable_file = 'GSM1647123_GM12878vsHL.indextable.txt'
    # readcounts_file = 'GSM1647123_GM12878vsHL.readcounts.txt'
    # indexconversion_file = 'GSM1647123_GM12878vsHL.indexconversion.txt'

    import yaml
    with open('file_list.yaml', 'r') as fopen:
        data = yaml.load(fopen)

    HL_path = data['main_dir']
    for sample in data['samples']:
        file_dir = list(sample.keys())[0]
        logger.info(file_dir)
        fastq_file = file_dir + '/' + sample[file_dir]['fastq_file']
        matrix_file = file_dir + '/' + sample[file_dir]['dhs_matrix']
        indexconversion_file = file_dir + '/' + sample[file_dir]['conversion']
        indextable_file = file_dir + '/' + sample[file_dir]['indextable']
        readcounts_file = file_dir + '/' + sample[file_dir]['readcounts']

        # add_tag_to_fastq(os.path.join(HL_path, fastq_file))
        logger.info('get barcodes templates list...')
        barcodes_templates_list = get_barcodes_templates(os.path.join(HL_path, indextable_file))

        logger.info('get fastq sequece counting dict...')
        barcodes_fastq_dict, barcodes_fastq_list = get_observed_barcodes(os.path.join(HL_path, fastq_file))

        logger.info('get index conversion table...')
        ic_list = get_index_conversion(os.path.join(HL_path, indexconversion_file))

        logger.info('get readcount barcodes...')
        rc_barcodes_list = get_readcount_barcodes(os.path.join(HL_path, readcounts_file))

        logger.info('get four different barcodes')
        split_barcodes_dict = split_barcodes(os.path.join(HL_path, indextable_file))

        logger.info('get distance between each barcodes')
        # get_ed_matrix(split_barcodes_dict)
        generate_indexconversion_table(split_barcodes_dict, barcodes_fastq_list, ic_list)
