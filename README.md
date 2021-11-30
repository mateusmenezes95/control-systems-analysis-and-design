# Control Systems Analysis and Design Scripts
Store scripts for the master's course Control Systems Analysis and Design

## Folder structure

The main folders and files are organized in the follow tree scheme:

```bash
├── doc
│   ├── first_report
│   │   ├── chapters
│   │   │   ├── challenge1
│   │   │   │   ├── challenge1.tex
│   │   │   │   ├── scenario1.tex
│   │   │   │   ├── scenario2.tex
│   │   │   │   ├── scenario3.tex
│   │   │   │   └── scenario4.tex
│   │   │   ├── challenge2
│   │   │   │   └── challenge2.tex
│   │   │   └── challenge3
│   │   │       └── challenge3.tex
│   │   ├── first_report.tex
│   │   └── images
│   │       ├── challenge1
│   │       │   ├── resultado-cenario1.tex
│   │       │   ├── resultado-cenario2-a.tex
│   │       │   ├── resultado-cenario2-b.tex
│   │       │   ├── resultado-cenario3.tex
│   │       │   ├── resultado-cenario4-a.tex
│   │       │   ├── resultado-cenario4-b.tex
│   │       │   ├── resultado-cenario4-c.tex
│   │       │   └── rlocus-question4.tex
│   │       ├── challenge2
│   │       │   ├── resultado-questao-1.tex
│   │       │   ├── resultado-questao-3-4.tex
│   │       │   ├── resultado-questao-5.tex
│   │       │   ├── resultado-questao-6-7.tex
│   │       │   ├── resultado-questao-6-9.tex
│   │       │   ├── resultado-questao-8.tex
│   │       │   └── resultado-questao-9.tex
│   │       ├── challenge3
│   │       │   ├── resultado-questao-2.tex
│   │       │   ├── resultado-questao-3-4-compensadores.tex
│   │       │   ├── resultado-questao-3-4-malha-aberta.tex
│   │       │   ├── resultado-questao-5-6-malha-aberta.tex
│   │       │   ├── resultado-questao-7-resposta-na-freq.tex
│   │       │   ├── resultado-questao-7-resposta-no-tempo.tex
│   │       │   ├── resultado-questao-8-9-compensadores.tex
│   │       │   ├── resultado-questao-8-9-malha-aberta.tex
│   │       │   └── resultado-questao-8-9-resposta-no-tempo.tex
│   │       └── challenge4
│   └── second_report
│       ├── chapters
│       │   ├── challenge4
│       │   │   ├── images
│       │   │   │   └── resultado-questao-5.tex
│       │   │   └── src
│       │   │       ├── challenge-4.tex
│       │   │       ├── imc.tex
│       │   │       ├── poles_allocation.tex
│       │   │       └── poles_allocation_vs_imc.tex
│       │   ├── challenge5
│       │   │   ├── images
│       │   │   └── src
│       │   │       └── challenge-5.tex
│       │   ├── challenge6
│       │   │   ├── images
│       │   │   │   ├── resultado-1-questao-10.tex
│       │   │   │   ├── resultado-2-questao-10.tex
│       │   │   │   ├── resultado-3-questao-10.tex
│       │   │   │   ├── resultado-4-questao-10.tex
│       │   │   │   ├── resultado-5-questao-10.tex
│       │   │   │   ├── resultado-6-questao-10.tex
│       │   │   │   ├── resultado-questao-7.tex
│       │   │   │   └── resultado-questao-9.tex
│       │   │   └── src
│       │   │       └── challenge-6.tex
│       │   ├── challenge7
│       │   │   ├── images
│       │   │   │   ├── resultado-1-questao-5.tex
│       │   │   │   ├── resultado-2-questao-5.tex
│       │   │   │   └── resultado-3-questao-5.tex
│       │   │   └── src
│       │   │       └── challenge-7.tex
│       │   └── initial_considerations.tex
│       └── second_report.tex
├── examples
│   ├── Simula_PMI_Matlab_Octave.m
│   └── Simula_PMI_sem_atraso_Octave.m
├── lib
│   ├── common_functions_script.m
│   └── graphs_functions_script.m
└── scripts
    ├── challenge-1
    │   ├── question-1.m
    │   ├── question-2.m
    │   ├── question-3.m
    │   └── question-4.m
    ├── challenge-2
    │   ├── common_parameters_scripts.m
    │   ├── question-1-4.m
    │   ├── question-5.m
    │   └── question-6-9.m
    ├── challenge-3
    │   ├── common_parameters_scripts.m
    │   ├── lead-compensator-4-6.m
    │   ├── question-2.m
    │   ├── question-3-4.m
    │   ├── question-5-6.m
    │   ├── question-7.m
    │   └── question-8-9.m
    ├── challenge-4
    │   ├── common_parameters_script.m
    │   ├── controllers
    │   ├── question-2-3.m
    │   ├── question-4.m
    │   └── question-5.m
    ├── challenge-6
    │   ├── base_script.m
    │   ├── question-6-10.m
    │   └── specific_state_space_graphs_funtions.m
    └── challenge-7
        └── question-1-5.m

```

## Dependencies

To run the scripts control and signal dependencies must be installed already. Whether not, inside octave-cli, run

```console
pkg install -forge control
pkg install -forge signal
pkg install -forge miscellaneous
```
## Tips

### Command to save graphs images in `.tex` extension

```console
print -depslatex -S"<width>,<height>" <image-path>/<image-name>
```
> Change names inside '<>' accordingly