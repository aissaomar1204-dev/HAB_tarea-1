# Tarea 1: Análisis Funcional 

## Descripción

Script de análisis de enriquecimiento funcional para los genes **COX4I2**, **ND1** y **ATP6** utilizando la API de g:Profiler. El análisis identifica procesos biológicos, funciones moleculares y componentes celulares significativamente enriquecidos.

## Resultados Principales

El análisis identificó **136 términos enriquecidos significativamente** distribuidos en:
- **65 procesos biológicos** (GO:BP)
- **23 funciones moleculares** (GO:MF)
- **21 componentes celulares** (GO:CC)
- **12 rutas KEGG**
- **9 rutas Reactome**
- **6 rutas WikiPathways**

Los genes están principalmente asociados con:
- Fosforilación oxidativa
- Cadena de transporte de electrones
- Respiración aeróbica
- Membrana interna mitocondrial

## Estructura del Proyecto

```
.
├── data/
│   └── genes_input.txt          # Genes de entrada (COX4I2, ND1, ATP6)
├── scripts/
│   └── analisis-funcional.py    # Script principal
├── results/
│   ├── resultados.tsv            # Resultados del análisis
│   └── TareaHAB_.pdf               # Informe detallado en LaTeX
├── requirements.txt              # Dependencias del proyecto
└── README.md                     # Este archivo
```

## Instalación y Ejecución

### Requisitos previos
- Python 3.7 o superior
- pip (gestor de paquetes de Python)

### Paso 1: Instalar dependencias

```bash
pip install -r requirements.txt
```

O en sistemas que requieran:

```bash
python -m pip install -r requirements.txt
```

### Paso 2: Ejecutar el análisis

Desde la carpeta `scripts/`:

```bash
python analisis_funcional.py -i ..\data\genes_input.txt -o ..\results\resultados.tsv
```

Desde la raíz del proyecto:

```bash
python scripts\analisis_funcional.py -i data\genes_input.txt -o results\resultados.tsv
```

## Salida

El script genera:

1. **resultados.tsv** - Archivo con los términos enriquecidos en formato TSV
2. **Salida en consola** - Resumen con los top 10 términos más significativos

### Formato del archivo TSV

| Columna | Descripción |
|---------|-------------|
| source | Base de datos (GO:BP, GO:MF, GO:CC, KEGG, REAC, WP) |
| native | ID del término |
| name | Nombre descriptivo del término |
| p_value | Valor p del test estadístico |
| significant | Significancia estadística (True/False) |
| intersection_size | Número de genes de entrada en el término |
| term_size | Tamaño total del término |
| query_size | Número total de genes analizados |
| intersections | Genes específicos presentes |

## Documentación Adicional

Para una explicación detallada del análisis, metodología e interpretación de resultados, consulte el documento **informe.pdf** en la carpeta `results/`.

## Metodología

- **Herramienta:** g:Profiler (https://biit.cs.ut.ee/gprofiler/)
- **Organismo:** Homo sapiens (hsapiens)
- **Test estadístico:** Test exacto de Fisher
- **Corrección:** FDR (False Discovery Rate)
- **Umbral de significancia:** p < 0.05

## Genes Analizados

| Gen | Nombre completo | Función |
|-----|----------------|---------|
| COX4I2 | Citocromo c oxidasa subunidad 4 isoforma 2 | Complejo IV de la cadena respiratoria |
| ND1 | NADH deshidrogenasa subunidad 1 | Complejo I de la cadena respiratoria |
| ATP6 | ATP sintasa F0 subunidad 6 | Complejo V (ATP sintasa) |

## Información Adicional

Para más información sobre g:Profiler: https://biit.cs.ut.ee/gprofiler/gost

---

**Autor:** Aissa Omar El Hammouti Chachoui  
**Fecha:** 4 de Noviembre 2025
