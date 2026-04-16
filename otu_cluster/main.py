import os
import bioepy as be
from multiprocessing import freeze_support


#######################################################################################
#################       AGGIORNALO PER IL TUO PERCORSO        #########################
FOLD_PATH = "C:/Users/andre/Desktop/my_project/combo_pf/prjct_gamma/" 
#FOLD_PATH = "C:/Users/admin/Desktop/my_project/combo_pf/prjct_gamma/"
#######################################################################################

KMER_LEN_MIN = 8
KMER_LEN_MAX = 12


if __name__ == "__main__":
    freeze_support()  # Supporto per Windows e Linux

 
    print ('This software implements a simplified reference-based OTU clustering algorithm for merged reads obtained from Paired End sequencing data. '
            'It follows these main steps:\n'
            'Merging PE-reads: The program merges paired-end (PE) sequencing reads (starting from two file.fastq) to '
            'reconstruct full-length sequences.\n'
            'Calculating Jaccard Distance: The merged sequences are then compared with a reference collection (reference_file.fasta). '
            'The Jaccard distance, a measure of similarity between sets, is used to quantify the differences '
            'between the merged sequences and reference sequences.\n'
            'Inferring Clusters: Based on the computed Jaccard distances, the software groups similar sequences '
            'into clusters. This helps in identifying related sequences and understanding the structure of the '
            'dataset.\n\n')
    
    #input ("Press enter to start...")

    # Scelta dei file per l'utente finale (versione per eseguibile)
    """
    FILE_PATH = input('Before to start the analysis put all your file you need in a fold and paste here its path: ').strip()
    while not os.path.isdir(FILE_PATH):
        FILE_PATH = input('\nThe directory selected does not exist. Check it and past here the correct one: ').strip()

    file1_name = input('\nPast here the name of the PE reads 1st file: ').strip()
    file1_path, file1_name = be.check_fastq_path(FOLD_PATH, file1_name)
    f1 = be.ParsingFastq(file1_path, file1_name)

    file2_name = input('\nPast here the name of the PE reads 2nd file: ').strip()
    file2_path, file2_name = be.check_fastq_path(FOLD_PATH, file2_name)
    f2 = be.ParsingFastq(file2_path, file2_name)

    fileref_name = input('\nPast here the name of the reference sequences file: ').strip()
    fileref_path, fileref_name = be.check_fastq_path(FOLD_PATH, file2_name)
    f2 = be.ParsingFasta(fileref_path, fileref_name)
    """

    file1_path = f'{FOLD_PATH}file/merging_test1.fq.gz'
    file2_path = f'{FOLD_PATH}file/merging_test2.fq.gz'
    fileref_path = f'{FOLD_PATH}file/ref.fa'
    f1 = be.ParsingFastq(file1_path, 'merging_test1.fq.gz')
    f2 = be.ParsingFastq(file2_path, 'merging_test2.fq.gz')

    fr = be.ParsingFasta(fileref_path, 'ref.fa')

    mpe = be.MergingPE(f1.valid_seqs, f2.valid_seqs)

    num_ref_seqs = len(fr.seqs)
    ref_seqs = [[fr.ids[i], fr.seqs[i]] for i in range (num_ref_seqs)]

    otu_cluster = be.ClusterOTU(FOLD_PATH, ref_seqs, mpe.merged_sequences)