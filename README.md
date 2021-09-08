# Buhlmann
This code is an implementation of the Buhlmann algorithm for decompression planning merged with IANTD standards. <br>
The decompression schedules predicted by this program are yet to be thoroughly tested in water. <br>
<br>
The author declines all responsability for any direct or indirect damage to objects, persons or animals deriving from the use of this program. <br>

## Authors
[Dr. Alberto Fabrizio](https://www.linkedin.com/in/alberto-fabrizio-03151b207/)

## Requirements

* `python >= 3.6`
* `numpy >= 1.16`
* `tabulate >= 0.8`

## General usage and options
Usage:
```
$python deco\_plan.py \[-h\] \[--altitude H\_ON\_SEA\] \[--alveolar ALVEOLAR\] --depth TDEPTH --time RUNT --fo2 FO2 \[FO2 ...\] --fhe FHE \[FHE ...\] \[--glow GF\_LOW\] \[--ghigh GF\_HI\] \[--last LAST\_DECO\] \[--debug\]
```

optional arguments:
  -h, --help           show this help message and exit <br>
  --altitude H\_ON\_SEA  Altitude on sea-level for the dive \[m\] \[default: 0.0\] <br>
  --alveolar ALVEOLAR  Which water vapor pressure in the lungs to use. \['buhl', 'schrein', 'navy'\] <br>
  --depth TDEPTH       Target maximum depth \[m\] <br>
  --time RUNT          Run time [min] at which you desire to quit the target maximum depth. <br>
  --fo2 FO2 [FO2 ...]  Fractions of oxygen in breathing gas (one value per tank). <br>
  --fhe FHE [FHE ...]  Fractions of helium in breathing gas (one value per tank). <br>
  --glow GF\_LOW        Gradient factor (Low) [default: 0.75]. <br>
  --ghigh GF\_HI        Gradient factor (High) [default: 0.75]. <br>
  --last LAST\_DECO     Last deco stop \[m\] \[default: 6\]. <br>
  --debug              Print Debug Info. <br>


## Examples of usage:

Dive at 40 m for 25 min bottom time carrying only air as a gas.

```
$ python deco_plan.py --depth 40 --time 25 --fo2 0.21 --fhe 0.0
```

Dive at 40 m for 25 min bottom time carrying air and EANx50.

```
$ python deco_plan.py --depth 40 --time 25 --fo2 0.21 0.5 --fhe 0.0 0.0
```

Dive at 40 m for 25 min bottom time carrying air and EANx50 with custom gradient factors.

```
$ python deco_plan.py --depth 40 --time 25 --fo2 0.21 0.5 --fhe 0.0 0.0 --glow 0.35 --ghigh 0.75
```

Dive at 40 m for 25 min bottom time carrying air and EANx50 with custom gradient factors at a altitude of 1000 m a.s.l.

```
$ python deco_plan.py --depth 40 --time 25 --fo2 0.21 0.5 --fhe 0.0 0.0 --glow 0.35 --ghigh 0.75 --altitude 1000
```

Dive at 50 m for 25 min bottom time carrying 20/45 trimix and EANx50

```
$ python deco_plan.py --depth 50 --time 25 --fo2 0.20 0.5 --fhe 0.45 0.0 
```

Dive at 50 m for 25 min bottom time carrying 20/45 trimix and EANx50 with last stop at 3 m


```
$ python deco_plan.py --depth 50 --time 25 --fo2 0.20 0.5 --fhe 0.45 0.0 --last 3
```








