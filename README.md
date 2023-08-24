# Detector-Esferico-Arrue

La parte mas importante es el archivo "construccion.cc". Ahí se encuentran todas las geometrias y la definición de los materiales con sus propiedades.
En el Construct() que se encuentra en este archivo se empieza definiendo todo el world y luego la semiesfera. Despues, depende del caso que se quiera estudiar, se construyen diferentes geometrias de detectores utilizando una sentencia de "if" con los messengers correspondientes. Si todos estan en false, se construye el detector mas simple. Despues de eso, se construyen los tubos que faltan y aparecen los distintos casos del teflon que tambien se activan con los messengers. Por último se construye la Mesh que por defecto tiene transparencia 1 (como si no estuviera), y se puede colocar el Xenon gaseoso en la tapa de la semiesfera opcionalmente.

En "generator.cc" se generá un fotón en una dirección aleatoria.

El archivo "dist.py" es el que se usa para crear las graficas apartir de los archivos que se crean con el macro "det.mac"."detgas.mac" es el macro equivalente cuando se coloca el Xenon gas y la distancia del PMT cambia.

