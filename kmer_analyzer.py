import collections
import math

def contar_km(seq, k):
    kmer_counts = collections.Counter()
    for i in range(len(seq) - k + 1):
        kmer = seq[i:i+k]
        kmer_counts[kmer] += 1
    return kmer_counts

""" def estadisticas(kmer_counts):
    total_kmers = sum(kmer_counts.values())
    unico = len(kmer_counts)
    
    if unico == 0:
        return total_kmers, unico, None, None

    mas_freq = kmer_counts.most_common(1)[0]
    menos_frec = kmer_counts.most_common()[-1] 
    
    return total_kmers, unico, mas_freq, menos_frec """

