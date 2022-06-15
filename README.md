En este trabajo, se realizaron los siguientes puntos:
* Metaanálisis de datos transcriptómicos de líneas de melanoma parentales tratadas con inhibidores de la vía MAPK
* Metaanálisis de datos transcriptómicos de líneas de melanoma resistentes a inhibidores de la vía MAPK resistentes 
* Generación de un clasificador que trate de distinguir células sensibles al tratamiento y resistentes.

A continuación, se describe cada uno de ellos.

## Metaanálisis de líneas tratadas

En este metaanálisis, se tomaron datos transcrptómicos de líneas de melanoma parentales tratadas con inhibidores de la vía MAPK durante 6, 8, 24 o 48 h. Los datos públicos de GEO utilizados fueron los siguientes:

- [GSE64741:](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE64741) Obenauf A375 6h MAPKi, Obenauf A375 48h MAPKi, Obenauf Colo800 48h y Obenauf UACC62 48h.
- [GSE87641:](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE87641) Fallahi Colo858 24h MAPKi, Fallahi Colo858 48h MAPKi, Fallahi MMACSF 24h MAPKi y Fallahi MMACSF 48h MAPKi.
- [GSE104869:](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE104869) Corre 501Mel 48h MAPKi.
- [GSE127988:](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE127988) Gerosa A375 24h MAPKi.
- [GSE141021:](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE141021) Smalley WM164 8h MAPKi.
- [GSE109731:](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE109731) Reganfendt A375 8h MAPKi.
- [GSE75313:](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE75313) Song M229 48h MAPKi y Song M238 48h MAPKi.

Estos datos de células parentales tratadas y sin tratar se descargaron de SRA y se obtuvieron las pseudocuentas para cada uno de los genes con salmon utilizando como referencia el [transcriptoma humano de Ensembl que deriva de la versión del genoma GRCh38/hg38](https://www.ensembl.org/Homo_sapiens/Info/Index). Los archivos de salida se cargaron en R y se procesaron. Los objetos procesados se pueden encontrar [aquí](https://drive.google.com/drive/folders/1ycnEiHrvURblxy61YuX_MlQDHroYPgNr?usp=sharing). Con estos datos, se realizaron análisis de expresión diferencial de tratamiento *vs.* control para cada línea y condición. Las tablas obtenidas se integraron mediante metaanálisis y se examinaron los resultados obtenidos. El código utlizado se encuentra en [treatment-MA](https://github.com/yberda/tfm-bioinfo/tree/main/treatment-MA). El archivo ["rankprod-treatment.R"](https://github.com/yberda/tfm-bioinfo/blob/main/treatment-MA/rankprod-treatment.R) contiene los metaanálisis por combinación de ránkings y ["treatment.R"](https://github.com/yberda/tfm-bioinfo/blob/main/treatment-MA/treatment.R), los metaanálisis por combinación de tamaño del efecto (TE-MA) y la  comparación y exploración de todos los resultados.

## Metaanálisis de líneas resistentes

En este metaanálisis, se tomaron datos transcriptómicos de líneas de melanoma resistentes a inhibidores de la vía MAPK. Los datos públicos de GEO utilizados fueron los siguientes:

- [GSE99923:](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE99923) Singleton A375R.
- [GSE108382:](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE108382) Ho 451LuR y Ho A375R.
- [GSE158119:](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE158119) Berico MM074R.
- [GSE75313:](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE75313) Song M238R, Song M229R Song M395R y Song SKMel28R.

Los pasos seguidos fueron equivalentes a los del apartado anterior. Los datos de salmon procesados que se utilizaron de partida se encuentran [aquí](https://drive.google.com/drive/folders/1ycnEiHrvURblxy61YuX_MlQDHroYPgNr?usp=sharing).

El código utilizado se encuentra en [resistance-MA](https://github.com/yberda/tfm-bioinfo/tree/main/resistance-MA). El archivo ["rankprod-resistance.R"](https://github.com/yberda/tfm-bioinfo/blob/main/resistance-MA/rankprod-resistance.R) contiene los metaanálisis por combinación de ránkings y ["resistance.R"](https://github.com/yberda/tfm-bioinfo/blob/main/resistance-MA/resistance.R), los metaanálisis por combinación de tamaño del efecto (TE-MA), el análisis de expresión diferencial de las líneas doble-resistentes y la  comparación y exploración de todos los resultados, incluida la comparación entre los metaanálisis de resistencia vs. los de tratamiento de parentales.

## Construcción del clasificador

Para la construcción del un clasificador capaz de distinguir células parentales y células resistentes, [se transformaron las cuentas de muestras de los dos casos en R](https://github.com/yberda/tfm-bioinfo/blob/main/classifier/tpm-norm.R) y se tomaron los valores obtenidos para el entrenamiento de varios modelos utilizando la biblioteca de `sci-kit learn` en Jupyter Notebook (python3). 

Los identificadores GEO de los datos utilizados son los siguientes: [GSE99923](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE99923), [GSE108382](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE108382), [GSE158119](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE158119), [GSE75313](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE75313), [GSE80829](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE80829), [GSE62526](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE62526), [GSE110054](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE110054), [GSE110948](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE110948), [GSE66539](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE66539), [GSE129127](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE129127), [GSE50535](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE50535), [GSE74729](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE74729), [GSE134459](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE134459) y [GSE65186](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE65186).

Los datos de salmon procesados que se utilizaron de partida se encuentran [aquí](https://drive.google.com/drive/folders/1ycnEiHrvURblxy61YuX_MlQDHroYPgNr?usp=sharing). El código utilizado para la transformación y el cuaderno ipynb con los clasificadores se encuentran en [classifier](https://github.com/yberda/tfm-bioinfo/tree/main/classifier). Tambén se puede ver el contenido del cuaderno a través del siguiente enlace: [cuaderno-clasificadores.](https://nbviewer.org/github/yberda/tfm-bioinfo/blob/main/classifier/classifier-cell-lines.ipynb)
