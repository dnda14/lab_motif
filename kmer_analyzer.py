import collections
import math

def contar_km(seq, k):
    kmer_counts = collections.Counter()
    for i in range(len(seq) - k + 1):
        kmer = seq[i:i+k]
        kmer_counts[kmer] += 1
    return kmer_counts

def estadisticas(kmer_counts):
    total_kmers = sum(kmer_counts.values())
    unico = len(kmer_counts)
    
    if unico == 0:
        return total_kmers, unico, None, None

    mas_freq = kmer_counts.most_common(1)[0]
    menos_frec = kmer_counts.most_common()[-1]
    
    return total_kmers, unico, mas_freq, menos_frec

def comparar(kmer_counts1, kmer_counts2):
    
    set1 = set(kmer_counts1.keys())
    set2 = set(kmer_counts2.keys())
    
    interseccion = set1.intersection(set2)
    union = set1.union(set2)
    
    jaccard_index = len(interseccion) / len(union) if len(union) > 0 else 0
    
    unique_to_seq1 = set1 - set2
    unique_to_seq2 = set2 - set1
    
    return jaccard_index, len(interseccion), len(unique_to_seq1), len(unique_to_seq2)
