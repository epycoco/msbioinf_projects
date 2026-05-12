import bioepy as be

FILE_FOLD_PATH = 'C:/Users/andre/Desktop/my_project/combo_pf/prjct_gamma/file'
#FILE_FOLD_PATH = 'C:/Users/admin/Desktop/my_project/combo_pf/prjct_gamma/file'


if __name__ == '__main__':
    chr_file_path = f'{FILE_FOLD_PATH}/chr7.fa'
    chr_file_name = 'chr7.fa'

    chr = be.ParsingFasta(chr_file_path, chr_file_name)

    while True:
        if chr.ids.__contains__('chromosome'):
            cftr = be.Transcript('cftr', 'ENSG00000001626', chr.ids, chr.seqs, f'{FILE_FOLD_PATH}/ensembl.txt.gz')
            break
        else:
            chr_file_path = input (f'The {chr_file_name} NOT contain a chromosome sequence.'
                                    'Check if you past it correctly or if the file is in the correct path '
                                    'and past here his name: ').strip()
            chr = be.ParsingFasta(chr_file_path, chr_file_name)




# ENSG00000001626



