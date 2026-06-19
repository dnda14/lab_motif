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
    """
    Selecciona la mejor posición evaluando k-mers individualmente.
    1. Para cada k-mer, calcula su mediana global.
    2. En cada secuencia, escoge la 1 posición más cercana a esa mediana.
    3. Calcula la dispersión (Desv. Est.) solo de esas posiciones filtradas.
    4. Elige el k-mer con menor dispersión como el verdadero motif.
    """
    import statistics
    
    print(f"\n{'='*60}")
    print(f"  ANÁLISIS DE DISPERSIÓN POR K-MER (1 Hit por Secuencia)")
    print(f"{'='*60}")
    print(f"  {'K-mer':<15} | {'Mediana':<10} | {'Dispersión (Stdev)':<20} | {'Seq_Hits'}")
    print(f"  {'-'*70}")

    estadisticas_kmers = [] # Guardará (kmer, dispersion, mejores_pos_por_seq, mediana)

    for kmer, _ in lista_candidatos:
        # 1. Encontrar TODAS las posiciones de este k-mer
        todas_pos = []
        pos_por_seq_raw = {}
        for k, v in set_secuencias.items():
            posiciones = []
            inicio = 0
            while True:
                pos = v.find(kmer, inicio)
                if pos == -1: break
                posiciones.append(pos)
                todas_pos.append(pos)
                inicio = pos + 1
            pos_por_seq_raw[k] = posiciones

        if not todas_pos:
            continue
            
        # Mediana inicial cruda para este k-mer
        mediana_cruda = statistics.median(todas_pos) #cambiar a mean si se requiere

        # 2. Elegir 1 sola posición por secuencia (la más cercana a la mediana cruda)
        posiciones_filtradas = []
        pos_elegida_por_seq = {}
        
        for k, posiciones in pos_por_seq_raw.items():
            if posiciones:
                mejor_pos = min(posiciones, key=lambda p: abs(p - mediana_cruda))
                posiciones_filtradas.append(mejor_pos)
                pos_elegida_por_seq[k] = mejor_pos
        
        # 3. Calcular la dispersión final de este k-mer con las posiciones filtradas
        n_hits = len(posiciones_filtradas)
        if n_hits > 1:
            mediana_final = statistics.median(posiciones_filtradas)# cambiar a mean si e requiere
            dispersion = statistics.stdev(posiciones_filtradas)
        elif n_hits == 1:
            mediana_final = posiciones_filtradas[0]
            dispersion = 0.0
        else:
            continue
            
        print(f"  {kmer:<15} | {mediana_final:<10.1f} | {dispersion:<20.2f} | n={n_hits}")
        estadisticas_kmers.append({
            'kmer': kmer,
            'dispersion': dispersion,
            'hits': n_hits,
            'mediana': mediana_final,
            'posiciones_seq': pos_elegida_por_seq
        })

    # 4. Elegir el k-mer ganador (prioridad: menor dispersión, luego más hits, luego más largo)
    # Ordenamos por: dispersión (ascendente), hits (descendente), longitud kmer (descendente)
    estadisticas_kmers.sort(key=lambda x: (x['dispersion'], -x['hits'], -len(x['kmer'])))
    
    ganador = estadisticas_kmers[0]
    kmer_ganador = ganador['kmer']
    mediana_ganadora = ganador['mediana']
    
    print(f"\n  >> K-MER GANADOR: {kmer_ganador} (Dispersión: {ganador['dispersion']:.2f}, Mediana: {mediana_ganadora:.1f})")
    
    # 5. Extraer regiones basadas en el k-mer ganador
    rangos_por_secuencia = {}
    for k, v in set_secuencias.items():
        if k in ganador['posiciones_seq']:
            mejor_pos = ganador['posiciones_seq'][k]
            
            # Límite más solidario: 20% de la longitud de la secuencia (ej. 160 pb -> 32 pb de margen)
            limite_distancia = int(0.20 * len(v))
            distancia = abs(mejor_pos - mediana_ganadora)
            
            if distancia <= limite_distancia:
                print(rellenar(kmer_ganador, v, mejor_pos) + f"  [Distancia: {distancia:.1f}]")
                rangos_por_secuencia[k] = (mejor_pos, mejor_pos + len(kmer_ganador))
            else:
                print(f"{k}: IGNORADA — La mejor coincidencia (dist {distancia:.1f}) supera margen del 20% ({limite_distancia} pb). Extrayendo a ciegas en el anclaje.")
                inicio_ciego = int(mediana_ganadora)
                rangos_por_secuencia[k] = (inicio_ciego, inicio_ciego + len(kmer_ganador))
        else:
            print(f"{k}: IGNORADA — No contiene el k-mer ganador. Extrayendo a ciegas en el anclaje.")
            inicio_ciego = int(mediana_ganadora)
            rangos_por_secuencia[k] = (inicio_ciego, inicio_ciego + len(kmer_ganador))
            
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
    k_values = [18] #18,19,25
    
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
    lista_candidatos = seleccionar_candidatos(k_values, len(secuencias_activas), umbral=1)
    
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
    # pregunta 6: matriz de frecuencias
    #--------
    
    print("\n" + "="*60)
    print("  6. MATRIZ DE FRECUENCIAS POR POSICIÓN")
    print("="*60)
    
    if not resultado_alineacion:
        print("No hay alineación para evaluar.")
        return

    longitud_alineamiento = len(resultado_alineacion[0])
    matriz_frecuencias = [Counter() for _ in range(longitud_alineamiento)]
    
    for seq in resultado_alineacion:
        for j in range(longitud_alineamiento):
            matriz_frecuencias[j][seq[j]] += 1
            
    # Imprimir matriz (A, C, G, T, -)
    bases = ['A', 'C', 'G', 'T', '-']
    header_bases = "".join([f"{b:>5}" for b in bases])
    print(f"{'Pos':>5} | {header_bases}")
    print("-" * 40)
    for j, conteos in enumerate(matriz_frecuencias):
        fila_conteos = "".join([f"{conteos.get(b, 0):>5}" for b in bases])
        print(f"{j:>5} | {fila_conteos}")
        
    #-------
    # pregunta 7 & 8: secuencia consenso y grado de conservación
    #--------
    print("\n" + "="*60)
    print("  7 & 8. SECUENCIA CONSENSO Y GRADO DE CONSERVACIÓN")
    print("="*60)
    
    consenso = ""
    posiciones_conservadas = []
    posiciones_variables = []
    total_secuencias = len(resultado_alineacion)
    
    for j, conteos in enumerate(matriz_frecuencias):
        bases_obs = [b for b, c in conteos.items() if c > 0]
        
        # Representación del consenso con la notación de corchetes del usuario
        if len(bases_obs) > 1:
            base_repr = f"[{''.join(sorted(bases_obs))}]"
        else:
            base_repr = bases_obs[0]
            
        consenso += base_repr
        
        # Seleccionar el nucleótido más frecuente para estadísticas
        base_mas_frecuente = conteos.most_common(1)[0]
        frecuencia_max = base_mas_frecuente[1]
        
        # Calcular porcentaje de conservación
        porcentaje = (frecuencia_max / total_secuencias) * 100
        
        if porcentaje == 100.0 and bases_obs[0] != '-':
            posiciones_conservadas.append(j)
            tipo = "Totalmente conservada (100%)"
        else:
            posiciones_variables.append((j, bases_obs))
            tipo = f"Variable {bases_obs} ({porcentaje:.1f}%)"
            
        print(f"  Pos {j:>2}: {base_repr:>5}  --> {tipo}")
        
    conservacion_global = (len(posiciones_conservadas) / longitud_alineamiento) * 100
    
    #-------
    # pregunta 9: Reporte final
    #--------
    print("\n" + "="*60)
    print("  9. REPORTE FINAL DEL MOTIF ENCONTRADO")
    print("="*60)
    print(f"  Longitud del motif: {longitud_alineamiento} posiciones")
    print(f"  Secuencia consenso: {consenso}")
    print(f"  Porcentaje global de conservación: {conservacion_global:.1f}%")
    print(f"  Número de secuencias que contienen el motif: {total_secuencias}")
    print("\n  Posición de aparición en cada secuencia:")
    for k, v in secuencias_activas.items():
        if k in rangos:
            i_seq, f_seq = rangos[k]
            print(f"    {k}: [{i_seq}:{f_seq}]")

if __name__ == "__main__":
    main()
