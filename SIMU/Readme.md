
# Molecular Dynamics in Python

### Launching code

It is probably easiest if you launch 2 terminal windows. The first to manage the jupyter notebooks
the second to run the python interpreter.

Start by examining the Jupyter notebook `Start.ipynb`
Navigate to its directory in a terminal and type

$ `jupyter-lab Start.ipynb`

In the second terminal window you can run the simulation

$ `python run.py`

### Format of reports

1. All binomes should have a working animation with N > 4 particles bouncing from walls and colliding pair-wise. This should run by executing `python run.py`.

2. You should also have attempted a project working as a binome (described in `Start.ipynb`). These should run from the command line, and have the name `project.py`. If you attempt more than one code create files `project1.py` etc. Every independent run should have its own file to launch it. We will not go into files to change options.

3. Each binome should prepare a report as a Jupyter-lab notebook. It should be limited to the equivalent of four printed pages in length. We favour concise and precise descriptions of what you did, the problems you met and the solutions that you favoured. Figures can be included by saving the figure within the Github directory, and using the syntax
```html
            <img src="name.png" width="400">   
```
to include the figure `name.png`.
This report should be called `Report.ipynb`. The notebook does NOT have to be dynamic. We will not execute it. French or English.

4. Mathematical notation can be used in Jupyter-lab, with the help of latex expressions.

5. The jupyter notebook should only be saved in the github account of One member of the binome. Write us a mail to say who is in the binome, and who saved the file to github.

6. If you have further instructions for us to properly run your codes, add them to the end of this file (Readme.md), which we will read on github.

### Notes from Students for teaching staff


- Pour exécuter la simulation avec les N petites particules (sans grosse particule) : exécuter dans le terminal le fichier run.py
- Pour exécuter la simulation avec les N petites particules ET la grosse particule : exécuter dans le termminal le fichier run_b.py
- La modification des paramètres rayons grosse/petites particules, masse grosse/petite particules, nombre de particules, taille de la boite pour le tracé des courbes dans la fonction plot se font directement dans le fichier plot.py : ligne 8 pour les courbes position/déplacement/déplacement quadratique en fonction du temps ; ligne 42 pour la simulation de distribtion des vitesses ; ligne 79 pour le déplacement quadratique en fonction du temps (nous n'avons pas redigé un seul et même set de paramètres  car chaque tracé à besoin de paramètres différents pour être exploité correctement)
- Pour tracer les différentes courbes présentent dans le rapport :
* Pour tracer la courbe abscisses/ordonnées de la grosse particule en fonction du temps : exécuter dans le terminal le fichier python -m src.plot.py en ayant précisé dans ce ficher "case = 1" et en étant dans le root du fichier.
* Pour tracer la courbe déplacement de la grosse particule en fonction du temps : exécuter dans le terminal le fichier python -m src.plot.py en ayant précisé dans ce ficher "case = 2" et en étant dans le root du fichier.
* Pour tracer la courbe déplacement quadratique moyen de la grosse particule en fonction du temps : exécuter dans le terminal le fichier python -m src.plot.py en ayant précisé dans ce ficher "case = 3" et en étant dans le root du fichier.
* Pour tracer la courbe déplacement quadratique moyen de la grosse particule en fonction du temps d'échantillonnage : exécuter dans le terminal le fichier python -m src.plot.py en ayant précisé dans ce ficher "case = 4" et en étant dans le root du fichier.
* Pour tracer l'histogramme de la distribution des vitesses de la grosse particule : exécuter dans le terminal le fichier python -m src.plot.py en ayant précisé dans ce ficher "case = 5" et en étant dans le root du fichier.
 
