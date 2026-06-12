import matplotlib.pyplot as plt
import os

def plot_top_kmers(kmer_conteo, top_n=20, output_archivo="top_kmers.png", titulo="Top K-mers"):
    if not kmer_conteo:
        return
    
    mas_comun = kmer_conteo.most_common(top_n)
    kmers, conteo = zip(*mas_comun)
    
    plt.figure(figsize=(10, 6))
    plt.bar(kmers, conteo, color='skyblue')
    plt.xlabel('K-mers')
    plt.ylabel('Frecuencia')
    plt.title(titulo)
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.savefig(output_archivo)
    plt.close()



