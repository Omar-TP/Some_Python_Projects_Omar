from src.brownian_motion import Brownian
from scipy.stats import norm
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import numpy as np

def position_big_part(index_begin,sample_time,repet): # Permet de récupérer les positions successives de la grosse particule en sortie de chaque itération de la fonction md_step, séparées à chaque fois par le temps d'échantillonage ; prend en entrées l'index de temps à partir duquel on récupère les positions de la grosse particule (pour pouvoir s'affranchir de la position cristalline des particules au départ), le temps d'échantillonnage et le nombre d'itération de la fonction md_step
    simulation = Brownian(sample_time, sigma=0.15,N=50,L=5,m=1,M=10,big_part_size=0.6) # Définition de la variable simulation, qui excécute la sous-classe Brownian
    time = [k*sample_time for k in range(repet)] # Création d'une liste de temps, multiple du temps d'échantillonage (chaque mesure de la position de la grosse particule est séparée par le temps d'échantillonage)
    posx = [] # Initialisation de la liste des abscisses des positions de la grosse particule
    posy = [] # Initialisation de la liste des abscisses des positions de la grosse particule
    for i in range(repet): # On itère "repet" fois la fonction md_step
        point =simulation.md_step()[0]
        posx.append(point[0]) # Récupération des abcisses de la grosse particule à la ième itération
        posy.append(point[1]) # Récupération des ordonnées de la grosse particule à la ième itération
    return(time[index_begin:],posx[index_begin:],posy[index_begin:]) # On retourne la liste de temps, abscisses et ordonnées à partir de l'index auquel on souhaite commencer à étudier la position de la grosse particule

def dep_big_part(index_begin,sample_time,repet): # Permet de calculer le déplacement de la grosse particule sur une durée égale au temps d'échantillonage ; mêmes paramètres d'entrée que la fonction précédente
    time,x,y = position_big_part(index_begin,sample_time,repet) # On récupère la liste de temps, abscisses et ordonnées obtenue avec la fonction précédente
    dep = [] # On initialise la liste des déplacements 
    for k in range(1,len(x)): # On commence à partir de 1 pour éviter les problèmes de "out of range"
        dep.append(np.sqrt((x[k]**2+y[k]**2))-np.sqrt((x[k-1]**2+y[k-1]**2))) # Calcul de la différence des normes des vecteurs position aux intants t et t + t_ech
    return(time,dep) # Renvoie la liste des temps et des déplacements

def dep_quad_big_part(index_begin,sample_time,repet): # Permet de calculer le déplacement quadratique de la grosse particule sur une durée égale au temps d'échantillonage ; mêmes paramètres d'entrée que la fonction précédente
    time,x,y = position_big_part(index_begin,sample_time,repet) # On récupère la liste de temps, abscisses et ordonnées obtenue avec la fonction position_big_part
    dep_quad = [] # On initialise la liste des déplacements quadratiques
    for k in range(1,len(x)):
        dep_quad.append(((np.sqrt((x[k]**2+y[k]**2)))-np.sqrt((x[k-1]**2+y[k-1]**2)))**2) # Calcul des déplacements quadratiques
    return(time,dep_quad) # Retourne la liste des temps et des déplacements quadratiques 

def dep_quad_moy_big_part(index_begin,sample_time,repet): # Permet de calculer le déplacement quadratique moyen de la grosse particule ; mêmes paramètres d'entrée que la fonction précédente
    time,x,y = position_big_part(index_begin,sample_time,repet) # On récupère la liste de temps, abscisses et ordonnées obtenue avec la fonction position_big_part
    dep_quad = [] # On initialise la liste des déplacements quadratiques
    dep_quad_moy = [] # On initialise la liste des déplacements quadratiques moyens
    for k in range(1,len(x)): # l'idée est de calculer le déplacement quadratique à un instant, puis de faire la moyenne sur tous les déplacements quadratiques aux instants qui précédent
        dep_quad.append((np.sqrt((x[k]**2+y[k]**2))-np.sqrt((x[k-1]**2+y[k-1]**2)))**2) # Calcul du déplacement quadratique et ajout à la liste dep_quad
        dep_quad_moy.append(np.mean(dep_quad)) # Moyenne à un instant t de la liste dep_quad --> déplacement quadratique moyen ajouté à la liste dep_quad_moy
    return(time,dep_quad_moy) # Retourne la liste des temps et des déplacements quadratiques moyens

def velocity_big_part(repet): # Permet de récupérer les vitesses successives de la grosse particule en sortie de chaque itération de la fonction md_step, séparées à chaque fois par le temps d'échantillonage ; prend en entrée de le nombre d'itération de la fonction md_step à réaliser 
    simulation = Brownian(sample_time=0.015, sigma=0.15,N=50,L=5,m=1,M=10,big_part_size=0.6) # Définition de la variable simulation, qui excécute la sous-classe Brownian
    vel_x = [] # Initialisation de la liste des vitesses selon x
    vel_y = [] # Initialisation de la liste des vitesses selon x
    for i in range(repet): 
        velocity =simulation.md_step()[1] # Récupération des vitesses lors de la ième itération
        vel_x.append(velocity[0]) # Récupération de la vitesse selon x de la grosse particule à la ième itération
        vel_y.append(velocity[1]) # Récupération de la vitesse selon y de la grosse particule à la ième itération 
    return(vel_x,vel_y) # Retourne les listes des vitesses selon x et y

index_begin = 0 # Définition de l'index de temps à partir duquel on commence à récupèrer les positions de la grosse particule
case = 5 # Définition du cas que l'on souhaite tracer : 1 = Abscisse et ordonnée de la grosse particule en fonction du temps ; 2 = Déplacement de la grosse particule en fonction du temps ; 3 = Déplacement quadratique moyen de la grosse particule en fonction du temps et fitting ; 4 = Déplacement quadratique moyen de la grosse particule en fonction du temps d'échantillonage et fitting ; 5 = Distribution des vitesses de la grosse particule et fitting Gaussien
if case == 1: # Position (abscisse et ordonnée en fonction du temps)
    time,x,y = position_big_part(index_begin,sample_time=0.015,repet=2000) # On récupère la liste de temps et les abscisses et ordonnées de la grosse particule 
    fig, (ax1, ax2) = plt.subplots(1, 2)
    ax1.plot(time,x) # Tracé de l'abscisse en fonction du temps
    ax1.set_xlabel('Temps')
    ax1.set_ylabel('x(t)')
    ax1.set_title("Évolution de l'abscisse du centre \n de la grosse particule \n en fonction du temps")
    ax1.grid()
    ax2.plot(time,y) # Tracé de l'ordonnée en fonction du temps
    ax2.set_xlabel('Temps')
    ax2.set_ylabel('y(t)')
    ax2.set_title("Évolution de l'ordonnée du centre \n de la grosse particule \n en fonction du temps")
    ax2.grid()
    plt.tight_layout()
    plt.show()

elif case == 2: # Déplacement en fonction du temps
    time,dep = dep_big_part(index_begin,sample_time=0.015,repet=2000) # On récupère la liste de temps et les déplacements de la grosse particule
    plt.plot(time[:-1],dep) # Tracé du déplacement en fonction du temps ; on enlève le dernier temps pour avoir le même nombre d'élement dans la liste de temps et dans la liste de déplacement
    plt.title('Déplacement de la grosse particule en fonction du temps')
    plt.xlabel('Temps')
    plt.ylabel(r'D(t,$\tau$)')
    plt.grid()
    plt.show()
    
elif case == 3: # Déplacement quadratique moyen en fonction du temps
    time,dep_quad_moy = dep_quad_moy_big_part(index_begin,sample_time=0.015,repet=2000) # On récupère la liste de temps et les déplacements quadratiques moyens de la grosse particule
    t_bis = np.linspace(time[0],time[-1],1000) # Création d'une liste de temps de 1000 points sur laquelle on va effectuer la régression linéaire
    reg = np.polyfit(time[:-1],dep_quad_moy,1) # Régression linéaire du déplacement quadratique en fonction du temps
    print(reg) # On affiche les coefficients de la régression linéaire
    reg_y = [reg[0]*t + reg[1] for t in t_bis] # Calcul des ordonnées de la droite de régression linéaire
    plt.plot(time[:-1],dep_quad_moy,'bx') # Tracé du déplacement quadratique moyen en fonction du temps ; on enlève le dernier temps pour avoir le même nombre d'élement dans la liste de temps et dans la liste de déplacement
    plt.plot(t_bis,reg_y,'b-') # Tracé de la régression linéaire
    plt.title('Déplacement quadratique moyen de la grosse particule \n en fonction du temps')
    plt.xlabel('Temps')
    plt.ylabel(r'$Q(t,\tau)$')
    plt.grid()
    plt.show()

elif case == 4: # Déplacement quadratique moyen en fonction du temps d'échantillonnage ; attention, cette portion met très longtemps à s'exécuter
    sample_time_list = [0.001, 0.00621053, 0.01142105, 0.01663158, 0.02184211, 0.02705263, 0.03226316, 0.03747368, 0.04268421, 0.04789474, 0.05310526, 0.05831579, 0.06352632, 0.06873684, 0.07394737, 0.07915789, 0.08436842, 0.08957895, 0.09478947, 0.1,1,4,7,10,20,50,100] # Définition de la liste des temps d'échantillonnages (on s'arrête à 100 car au dessus, cela prend beaucoup trop de temps)
    moyenne = [] # Initialisation de la liste des déplacements quadratiques moyens
    for sample in sample_time_list: # On parcours la liste des temps d'échantillonages
        time,dep_quad = dep_quad_big_part(index_begin,sample_time=sample,repet=10) # Calcul du déplacement quadratique moyen
        moyenne.append(abs(np.mean(dep_quad))) # Calcul de la moyenne de la liste des déplacements quadratiques moyens au bout de 10 itérations de md_step
    plt.xscale("log") # Échelle log en abscisse
    plt.yscale("log") # Échelle log en ordonnée
    plt.plot(sample_time_list,moyenne,'ro') # Tracé du déplacement quadratique moyen en fonction du temps d'échantillonage
    plt.title("Déplacement quadratique moyen en fonction du temps d'échantillonage")
    plt.xlabel(r"Temps d'échantillonage $\tau$")
    plt.ylabel('Déplacement quadratique moyen')
    plt.grid()
    plt.show()

elif case == 5: # Distribution des vitesses de la grosse particule
    vel_x,vel_y = velocity_big_part(100) # On récupère les vitesses selon x et selon y de la grosse particule
    (mux, sigmax) = norm.fit(vel_x) # Fittin Gaussien de la vitesse selon x, et calcul des moyennes et écarts types de la Gaussienne fittée
    (muy, sigmay) = norm.fit(vel_y) # Fitting Gaussien de la vitesse selon y, et calcul des moyennes et écarts types de la Gaussienne fittée
    fig, (ax1, ax2) = plt.subplots(1, 2)
    nx, binsx, patchesx = ax1.hist(vel_x, 50, density='True') # Tracé de l'histogramme des vitesses selon x avec 50 rectangles
    ny, binsy, patchesy = ax2.hist(vel_y, 50, density='True') # Tracé de l'histogramme des vitesses selon y avec 50 rectangles
    yx = norm.pdf( binsx, mux, sigmax) # Calcul des ordonnées de la régression Gaussienne de la distribution des vitesses selon x
    yy = norm.pdf( binsy, muy, sigmay) # Calcul des ordonnées de la régression Gaussienne de la distribution des vitesses selon y
    lx = ax1.plot(binsx, yx, 'r--', linewidth=2) # Tracé de la régression Gaussienne pour les vitesses selon x
    ax1.set_title('Distribution de la vitesse \n selon x normalisée')
    ax1.set_xlabel('Vitesse selon x')
    ax1.set_ylabel('Densité')
    ly = ax2.plot(binsy, yy, 'r--', linewidth=2) # Tracé de la régression Gaussienne pour les vitesses selon y
    ax2.set_title('Distribution de la vitesse \n selon y normalisée')
    ax2.set_xlabel('Vitesse selon y')
    ax2.set_ylabel('Densité')
    plt.tight_layout()
    plt.show()
    print('Selon x', 'mu =', mux,'sigma =',sigmax) # Affichage des écart type et moyenne de la régression Gaussienne de vx
    print('Selon y','mu =',muy,'sigma =',sigmay) # Affichage des écart type et moyenne de la régression Gaussienne de vy



