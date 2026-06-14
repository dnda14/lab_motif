import kmer_analyzer
import visualizador
from collections import Counter,defaultdict
import alineacion_multiple

set_frecuencias_seq_globales=defaultdict(Counter)

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
            

def rellenar(subsecuencia,secuencia,pos):
    lon = len(secuencia)
    lon_sub = len(subsecuencia)
    prefijo = "."*pos
    sufijo = "."*(lon-(lon_sub+pos))
    
    return prefijo +subsecuencia +sufijo
    
def determinar_posiciones(set_secuencias,lista_candidatos):
    i_region = 25
    f_region = -1
    for k,v in set_secuencias.items():
        print(k+":")
        print(v)
        for i in lista_candidatos:
            pos = v.find(i[0])
            if pos>=0:
                i_region=min(i_region,pos)
                f_region=max(f_region,pos+len(i[0]))
                print(rellenar(i[0],v,pos))
    
    return i_region,f_region
            
            
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
    
    #---------
    # pregunta 2 : kmers candidatos
    #---------
    for k,v in set_secuencias.items():
        analiz_seq(v, k_values, k)
    
    #print(set_frecuencias_seq_globales)
    
    visualizador.plot_top_kmers(set_frecuencias_seq_globales[7], top_n=20,output_archivo=f"histogramas/{7}_top_frec_k{7}.png" , titulo="Top 7-mers")
    visualizador.plot_top_kmers(set_frecuencias_seq_globales[8], top_n=20,output_archivo=f"histogramas/{8}_top_frec_k{8}.png" , titulo="Top 8-mers")
    visualizador.plot_top_kmers(set_frecuencias_seq_globales[9], top_n=20,output_archivo=f"histogramas/{9}_top_frec_k{9}.png" , titulo="Top 9-mers")
    
    #visualmente se ven 6 candidatos(3 en 7mer, 2 en 8mer y 1 en 9mer) todos con presencia en las 10 secuencias, las otros posibles candidatos tienen solo  prsencia en 4 secuencias en todos los 3 kmern disponibles.
    lista_candidatos = []
    
    lista_candidatos+=set_frecuencias_seq_globales[7].most_common(4)
    lista_candidatos+=set_frecuencias_seq_globales[8].most_common(3)
    lista_candidatos+=set_frecuencias_seq_globales[9].most_common(2)
    
    #print(lista_candidatos)
    
    #-------
    # pregunta 3: localizar ocurrencias
    #--------
    i_region,f_region = determinar_posiciones(set_secuencias,lista_candidatos)

    
    #-------
    # pregunta 4: extraer region conservada 
    #--------
    # par ala extraccion de la region se considera un rango maximo que tengan en comun todas las secuencias, ene este caso un rango con el k = 9, como es es unico se tomara ese kmer como la regiona conservada
    
    #tomaremos una region en comun donde se contenga las secuencias en comun entre los 10 y entre 4, porque ahi se ve un mayor peso de freceuncias de subsecuecnias similares entre las 10 seceuncias.
    
    set_regiones =[]
    print(f"la region abarca entre el indice{i_region} y {f_region}")
    
    for k,v in set_secuencias.items():
        set_regiones.append((k+'_r',v[i_region:f_region]))
        print(k+" :",end='\t')
        print(v[i_region:f_region])
        
    #-------
    # pregunta 5:alineamiento multiple 
    #--------
    
    resultado_alineacion =alineacion_multiple.run_alineacion(set_regiones)
    
    #-------
    # pregunta 6: matriz de recuecnias
    #--------
    
    matriz_frec = {}
    for i in resultado_alineacion:
        for j in range(len(i)):
            if j in matriz_frec:
                matriz_frec[j].add(i[j])
            else:
                matriz_frec[j]={i[j]}
    #print(matriz_frec)
    
    for k,v in matriz_frec.items():
        print(f"{k}\t|  {v}")
    
    #----------
    # pregunta  7: secuencia consenso
    #--------
    consenso=""
    for v in matriz_frec.values():
        consenso+=''.join(['[']+list(v)+[']']if len(list(v))>1 else v)
    print(consenso)
        
    
    
if __name__ == "__main__":
    main()
