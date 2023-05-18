# Photometry reduction
---

Данная программа предназначена для перевода фотометрических изображений к стандартным фотометричесским системам. На данный момент доступны стандартные каталоги: NOMAD, SDSS, PanSTARRS.
Программа способна работаь в 2х основных режимах:
- Изображение в одном фильтре. В этом случае происходит линейное преобразование к стандартной системе, без корекции цветов.
- Обработка изображения в 3х цветовых полосах, что позволяет произвести цветовую коррекцию изображений.

Апертураная фотометрия проводилась при помощи пакета [photutils](https://photutils.readthedocs.io/en/stable)

## Установка:
Программа написана на языке [python 3](https://www.python.org).

Также для корректной работы необходимо установить [ds9](https://sites.google.com/cfa.harvard.edu/saoimageds9/home?authuser=0).

Ниже приведены используемы python пакеты:

```
pip3 install photutils astropy pathlib numpy scipy matplotlib
```

## Пример использования программы:

Ниже приведена последовательность пунктов, которые рекомендовано учесть при запуске программы:

1. Программа запускается через основной скрипт main.py

 `python3 main.py`

2. В качестве обязательных аргументов необходимо подать на вход изображения:

`python3 main.py first_image.fits second_image.fits third_image.fits` --- для 3х цветового режима 

или

`python3 main.py first_image.fits N N` --- для монохроматического режима.
 

3. Далее посредством ключа `-с` выбирается один из каталогов: _NOMAD, SDSS, PS1_. (По умолчанию NOMAD)

`python3 main.py first_image.fits second_image.fits third_image.fits -c NOMAD`

4. Теперь необходимо выбрать цветовые полосы для каждого изображения. Указывать их следует слитно, без пробелов. Например, BVR. ОБРАТИТЕ ВНИМАНИЕ, что при запуске монохромного режима используется ПЕРВЫЙ из указываемых фильтров. Все последующие символы будут проигнорированы. (По умолчанию BVR)

`python3 main.py first_image.fits second_image.fits third_image.fits -c NOMAD -f BVR`

или 

`python3 main.py N N -c NOMAD -f R`

5. Для проведения корректной фотометрии стоит брать яркие, но не перенасыщенные изображения звёзд. Для регулировки есть два ключа `-lm` и `-um`. Они задают диапазон $LM \le m_{obs} \le UM$ из которого будут взяты стандарты. По умолчанию диапазон [14; 17]

`python3 main.py first_image.fits second_image.fits third_image.fits -c NOMAD -f BVR -lm 14 -um 17`

6. Ключом `--fwhm` задаётся оценочная FWHM звезды в пикселях. Её, к примеру, можно измерить при помощи программы ds9. Данная величина необходима для задания апретур. В данной программе используются круглые апертуры. Их радиус регулируется ключом `--k0`: $r_{aper} = k0*FWHM (pix)$. Для локальной оценки фона используются круглые кольца, радиусами $r_[in} = k1*FWHM (pix)$ и $r_{out}=k2*FWHM (pix)$. Они аналогично регулируются ключами `--k1` и `--k2`. (По умолчанию FWHM = 9 pix, k0 = 2; k1 = 2.5; k2 = 3)

`python3 main.py first_image.fits second_image.fits third_image.fits -c NOMAD -f BVR -lm=14 -um=17 --fwhm=9 --k0=2 --k1=2.5 --k2=3`

7. Также в программе доступен ручной режим посредством вызова ds9, позволяюший при необходимости вручную регулировать размеры апертур, а также удалять апертуры плохих источников. Подробнее о ручном режиме будет изложено ниже. Активация ручного режима происхожит при помощи ключа `--manual`

`python3 main.py first_image.fits second_image.fits third_image.fits -c NOMAD -f BVR -lm=14 -um=17 --fwhm=9 --k0=2 --k1=2.5 --k2=3 --manual`

8. ПО умолчанию все результаты работы программы сохраняются в той же дирректории, в которой был запущен скрипт. Для смены пути используйте ключ `-out`: 

`python3 main.py first_image.fits second_image.fits third_image.fits -c NOMAD -f BVR -lm=14 -um=17 --fwhm=9 --k0=2 --k1=2.5 --k2=3 --manual -out=path_to_output_dirrectory`

9. При обработке изображений порой возникает необходимость отбросить объекты, близкие к краю изображения. Для этого существует параметр `-bdr`, задающий в пикселях рамку по краям всего изображения, данные из которой не будут учитываться при обработке. Отмечу, что объект считается внутри рамки, если координаты его центра лежат в ней. (По умолчанию bdr=100)

`python3 main.py first_image.fits second_image.fits third_image.fits -c NOMAD -f BVR -lm=14 -um=17 --fwhm=9 --k0=2 --k1=2.5 --k2=3 --manual -out=path_to_output_dirrectory -bdr=100`

10. Если возникает необходимость замаскировать часть изображения, то для этого существует ключ `-p2m` В нем узазывется путь до FITS  файла маски.

`python3 main.py first_image.fits second_image.fits third_image.fits -c NOMAD -f BVR -lm=14 -um=17 --fwhm=9 --k0=2 --k1=2.5 --k2=3 --manual -out=path_to_output_dirrectory -bdr=100 -p2m mask.fits`

Мы описали все ключи, используемые при запуске данной программы. Для получения информации во время работы используйте ключ `-h`

## Результат работы программы

В результате работы программы в итоге имеются:
- FITS файлы изображений, приведенные к стандартной фотометрической системе.
- Рисунок, отображающий звёзды, исполбзованные для колибровки
- Регионовские файлы DS9, используемых апертур
- Таблица с результатами фотометрии и стандартами
- Уравнения перехода, полученные по обработанным данным.

## Алгоритм работы программы
![Block diagram with the program algorithm](diagram.png)
