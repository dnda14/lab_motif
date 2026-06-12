import os
import time
import tracemalloc
import kmer_analyzer
import visualizador
from collections import Counter,defaultdict

set_frecuencias_seq_globales=defaultdict(Counter)
results = []
def analiz_seq(seq, k_values, prefix):
    
    if not seq:
        print("No se pudo leer la secuencia.")
        return None
    
    print(f"Longitud de la secuencia: {len(seq)} pb")
    
    for k in k_values:
        print(f"> Extrayendo y contando {k}-mers...")
        
        kmer_counts = kmer_analyzer.contar_km(seq, k)
                
        #print(f"{k}-mers únicos: {len(kmer_counts)}")
        
        set_frecuencias_seq_globales[k] += kmer_counts
                
        visualizador.plot_top_kmers(kmer_counts, top_n=20,output_archivo=f"histogramas/{prefix}_histograma_k{k}.png" , titulo="Top K-mers")
            
    return results

def main():
    
    set_secuencias={"S1" : "ATCGTACGATGACCTGATCG",
                    "S2" : "GGTATACGATGACGTTACCA",
                    "S3" : "TTTCTACGATGACCATAGGT",
                    "S4" : "AACGTACGATGACGGGTTAA",
                    "S5" : "CGGATACGATGACTTCCGTA",
                    "S6" : "TACCTACGATGACAGGTACA",
                    "S7" : "GACTTACGATGACCGATAGC",
                    "S8" : "TCGATACGATGACTGGCAAT",
                    "S9" : "AGGCTACGATGACATTCGGA",
                    "S10" : "CCTATACGATGACGGAATTC"}
    
    
    k_values = [7, 8, 9]
    
    lista_resultados = []
    for k,v in set_secuencias.items():
        lista_resultados.append(analiz_seq(v, k_values, k))
    
    #print(set_frecuencias_seq_globales)
    
    visualizador.plot_top_kmers(set_frecuencias_seq_globales[7], top_n=20,output_archivo=f"histogramas/{7}_top_frec_k{7}.png" , titulo="Top 7-mers")
    visualizador.plot_top_kmers(set_frecuencias_seq_globales[8], top_n=20,output_archivo=f"histogramas/{8}_top_frec_k{8}.png" , titulo="Top 8-mers")
    visualizador.plot_top_kmers(set_frecuencias_seq_globales[9], top_n=20,output_archivo=f"histogramas/{9}_top_frec_k{9}.png" , titulo="Top 9-mers")
    
    set_frecuencias_seq_globales[7].most_common()
    set_frecuencias_seq_globales[8].most_common()
    set_frecuencias_seq_globales[9].most_common()

if __name__ == "__main__":
    main()
