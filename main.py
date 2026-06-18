import kmer_analyzer
import visualizador
from collections import Counter,defaultdict
import alineacion_multiple

# Frecuencia global (suma de conteos) por k
set_frecuencias_seq_globales=defaultdict(Counter)
# Presencia: para cada k, guarda un dict {kmer: set de nombres de secuencias donde aparece}
presencia_por_secuencia = defaultdict(lambda: defaultdict(set))

def analiz_seq(seq, k_values, prefix):
    
    if not seq:
        print("No se pudo leer la secuencia.")
        return None
    
    print(f"Longitud de la secuencia: {len(seq)} pb")
    
    for k in k_values:
        print(f"> Extrayendo y contando {k}-mers...")
        kmer_counts = kmer_analyzer.contar_km(seq, k)
        
        set_frecuencias_seq_globales[k] += kmer_counts
        
        # Registrar en cuáles secuencias aparece cada kmer
        for kmer in kmer_counts:
            presencia_por_secuencia[k][kmer].add(prefix)
                
        visualizador.plot_top_kmers(kmer_counts, top_n=20,output_archivo=f"histogramas/{prefix}_histograma_k{k}.png" , titulo="Top K-mers")


def seleccionar_candidatos(k_values, total_secuencias, umbral=0.3):
    """
    Selecciona k-mers candidatos de forma automática.
    
    Criterio: un k-mer es candidato si aparece en al menos
    (umbral * total_secuencias) secuencias distintas.
    
    Luego elimina k-mers redundantes (substrings de otro candidato más largo).
    
    Retorna lista de tuplas (kmer, n_secuencias_presentes) ordenada por
    longitud descendente y luego por presencia descendente.
    """
    min_presencia = int(umbral * total_secuencias)
    
    # 1. Recolectar todos los k-mers que superen el umbral de presencia
    todos_candidatos = {} # {kmer: n_secuencias}
    for k in k_values:
        for kmer, seqs in presencia_por_secuencia[k].items():
            n_seqs = len(seqs)
            if n_seqs >= min_presencia:
                todos_candidatos[kmer] = n_seqs
    
    print(f"\n{'='*60}")
    print(f"  SELECCIÓN AUTOMÁTICA DE CANDIDATOS")
    print(f"  Umbral: presencia en >= {min_presencia}/{total_secuencias} secuencias ({umbral*100:.0f}%)")
    print(f"{'='*60}")
    print(f"  K-mers que superan el umbral: {len(todos_candidatos)}")
    for kmer, n in sorted(todos_candidatos.items(), key=lambda x: (-len(x[0]), -x[1])):
        print(f"    {kmer} (len={len(kmer)}, presente en {n}/{total_secuencias} seqs)")
    
    # 2. Eliminar k-mers redundantes (substrings de otro candidato más largo)
    kmers_ordenados = sorted(todos_candidatos.keys(), key=len, reverse=True)
    no_redundantes = []
    for kmer in kmers_ordenados:
        es_substring = False
        for ya_incluido in no_redundantes:
            if kmer in ya_incluido:
                es_substring = True
                break
        if not es_substring:
            no_redundantes.append(kmer)
    
    # 3. Construir resultado final
    resultado = [(kmer, todos_candidatos[kmer]) for kmer in no_redundantes]
    resultado.sort(key=lambda x: (-len(x[0]), -x[1]))
    
    print(f"\n  Candidatos finales (sin redundancia):")
    for kmer, n in resultado:
        print(f"    {kmer} (len={len(kmer)}, presente en {n}/{total_secuencias} seqs)")
    
    return resultado
            

def rellenar(subsecuencia,secuencia,pos):
    lon = len(secuencia)
    lon_sub = len(subsecuencia)
    prefijo = "."*pos
    sufijo = "."*(lon-(lon_sub+pos))
    
    return prefijo +subsecuencia +sufijo
    
def determinar_posiciones(set_secuencias,lista_candidatos):
    rangos_por_secuencia = {}
    for k,v in set_secuencias.items():
        print(k+":")
        print(v)
        i_seq = len(v)
        f_seq = -1
        for i in lista_candidatos:
            pos = v.find(i[0])
            if pos>=0:
                i_seq=min(i_seq,pos)
                f_seq=max(f_seq,pos+len(i[0]))
                print(rellenar(i[0],v,pos))
        
        if f_seq != -1:
            rangos_por_secuencia[k] = (i_seq, f_seq)
        else:
            rangos_por_secuencia[k] = (0, len(v))
    
    return rangos_por_secuencia
            
            
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
    set_secuencias2 = {
                    "S1" : "ATCGGCTAACGTAGCTAGCTTGACCGTACGATCGATCGGATCGTAGCTAGCATCGATCGTACGATCGATGCTAGCTAGCATCGATCGATACGATCGTAGCTAGCTACGTAGCTAGCTACGTAGCTTACGATGACGGTACCGATCGATCGTAGCTAACGTA",
                    "S2" : "GCTAGCTAGCATCGATCGTAGCTAGCTAGCGATCGTAGCTAGCATCGATCGATCGTAGCTAGCTAGCATCGATCGATCGTAGCTAGCTAGCATCGATACGTAGCTACGTACGATGACATCGTAGCTAGCTAACGTAGCTAGCTAGCGATCGTAGCTAGCTA",
                    "S3" : "CGTAGCTAGCTAGCATCGATCGTAGCTAACGTAGCTAGCATCGATCGTAGCTAGCTAGCATCGATCGATCGTAGCTAGCTAACGTAGCTAGCTAGCATCGATCGTAGCTAGCTACGTACGATGATCGATCGTAGCTAACGTAGCTAGCTAGCATCGATCGTA",
                    "S4" : "AACGTAGCTAGCTAGCATCGATCGTAGCTAGCTAACGTAGCTAGCTAGCATCGATCGTAGCTAGCATCGATCGTAGCTAACGTAGCTAGCTAGCATCGATCGTAGCTAGCTACGTATACGATGACGCTAGCTAACGTAGCTAGCATCGATCGTAGCTAACGT",
                    "S5" : "TAGCTAGCATCGATCGTAGCTAACGTAGCTAGCTAGCATCGATCGTAGCTAACGTAGCTAGCATCGATCGTAGCTAACGTAGCTAGCTAGCATCGATCGTAGCTAGCTAACGTATACGATGTCCGTAGCTAACGTAGCTAGCTAGCATCGATCGTAGCTAACG",
                    "S6" : "GATCGTAGCTAACGTAGCTAGCATCGATCGTAGCTAACGTAGCTAGCTAGCATCGATCGTAGCTAACGTAGCTAGCATCGATCGTAGCTAACGTAGCTAGCTAGCATCGATCGTAGCTATACGATGACGATCGTAGCTAACGTAGCTAGCATCGATCGTAGCT",
                    "S7" : "CTAGCTAACGTAGCTAGCATCGATCGTAGCTAACGTAGCTAGCATCGATCGTAGCTAACGTAGCTAGCTAGCATCGATCGTAGCTAACGTAGCTAGCATCGATCGTAGCTAACGTATACGATGCCGCTAACGTAGCTAGCATCGATCGTAGCTAACGTAGCTA",
                    "S8" : "CGATCGTAGCTAACGTAGCTAGCATCGATCGTAGCTAACGTAGCTAGCTAGCATCGATCGTAGCTAACGTAGCTAGCATCGATCGTAGCTAACGTAGCTAGCTAGCATCGATCGTAGCTATACGATGACCGTAGCTAACGTAGCTAGCATCGATCGTAGCTAAC"
                    
    }
    k_values = [18,19,25]
    
    secuencias_activas = set_secuencias2  
    
    #---------
    # pregunta 2 : kmers candidatos
    #---------
    for k,v in secuencias_activas.items():
        analiz_seq(v, k_values, k)
    
    #print(set_frecuencias_seq_globales)
    
    for k in k_values:
        # Plot 1: frecuencia bruta (total de ocurrencias sumadas)
        visualizador.plot_top_kmers(set_frecuencias_seq_globales[k], top_n=20,output_archivo=f"histogramas/{k}_top_frec_k{k}.png" , titulo=f"Top {k}-mers (frecuencia total)")
        
        # Plot 2: presencia en secuencias distintas
        presencia_counter = Counter({kmer: len(seqs) for kmer, seqs in presencia_por_secuencia[k].items()})
        visualizador.plot_top_kmers(presencia_counter, top_n=20, output_archivo=f"histogramas/{k}_presencia_k{k}.png", titulo=f"Top {k}-mers (presencia en secuencias)")
    
    # Selección automática: filtra por presencia en secuencias distintas
    # y elimina k-mers redundantes (substrings de uno más largo)
    lista_candidatos = seleccionar_candidatos(k_values, len(secuencias_activas), umbral=0.8)
    
    #-------
    # pregunta 3: localizar ocurrencias
    #--------
    rangos = determinar_posiciones(secuencias_activas,lista_candidatos)

    
    #-------
    # pregunta 4: extraer region conservada 
    #--------
    # par ala extraccion de la region se considera un rango maximo que tengan en comun todas las secuencias, ene este caso un rango con el k = 9, como es es unico se tomara ese kmer como la regiona conservada
    
    #tomaremos una region en comun donde se contenga las secuencias en comun entre los 10 y entre 4, porque ahi se ve un mayor peso de freceuncias de subsecuecnias similares entre las 10 seceuncias.
    
    set_regiones =[]
    print("Rangos extraidos por secuencia:")
    
    for k,v in secuencias_activas.items():
        i_seq, f_seq = rangos[k]
        region_extraida = v[i_seq:f_seq]
        set_regiones.append((k+'_r',region_extraida))
        print(f"{k} [{i_seq}:{f_seq}] :\t{region_extraida}")
        
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
