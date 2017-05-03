'''
kmer_enrichment.py

counts the occurrences of all sequences length k in a given set of sequences.
If given a second file, will calculate the enrichment of File_A to File_B
for each k-mer.
** does NOT return a count for each enhancer

kmer_enrichment.py PATH/TO/FASTA_A k (PATH/TO/FASTA_B)
'''
import sys
import string
import itertools

def main():
    print '#',' '.join(sys.argv)
    f = [line.strip() for line in open(sys.argv[1], 'r')]
    k = int(sys.argv[2])
    try:
        enr = True
        f_2 = [line.strip() for line in open(sys.argv[3], 'r')]
    except:
        enr = False

#   populate kmers with keys
    kmers = {}
    kmers_2 = {}
    for kmer in itertools.product('ACGT',repeat=k):
        kmers[''.join(kmer)] = 0
        if enr:
            kmers_2[''.join(kmer)] = 0

#   count occurrences in file a
    seq = ""
    n = -1 # because it picks up an extra counting the last sequence in the file
    for i in range(len(f)):
        if f[i].startswith('#'): continue
        if f[i].startswith('>') or i == len(f)-1:
            n += 1
            bp_ind = 0
            while bp_ind < (len(seq)-(k-1)):
                kmers[string.upper(seq[bp_ind:bp_ind+k])] += 1
                bp_ind += 1
            seq = ""
            continue
        seq += f[i]

# calculate and print enrichment
    if enr:
        seq = ""
        n_2 = -1 # because it picks up an extra counting the last sequence in the file
        for i in range(len(f_2)):
            if f_2[i].startswith('#'): continue
            if f_2[i].startswith('>') or i == len(f_2)-1:
                n_2 += 1
                bp_ind = 0
                while bp_ind < (len(seq)-(k-1)):
                    kmers_2[string.upper(seq[bp_ind:bp_ind+k])] += 1
                    bp_ind += 1
                seq = ""
                continue
            seq += f_2[i]
        print "# k = %d; N_1 = %d vs N_2 = %d\nkmer\tenr" % (k,n,n_2)
        for key in sorted(kmers.keys()):
            # e = mean(FILE_A)/mean(FILE_B)
            try: e = (float(kmers[key])/n)/(float(kmers_2[key])/n_2)
            except: e = -1.0
            print "%s\t%.3f" % (key,e)
# if enrichment wasn't a thing, just print the counts for File_A
    else:
        print "# k = %d; N = %d\nkmer\tcount" % (k,n)
        for key in sorted(kmers.keys()):
            print "%s\t%d" % (key,kmers[key])

if __name__ == "__main__":
    main()
