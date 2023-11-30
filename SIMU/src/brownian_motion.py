from src.simul import Simul
import numpy as np
import matplotlib.pyplot as plt

def remove_center_particles(config_pos,config_vel,config_sigma,config_masse,center_pos,big_part_size): # On souhaite pouvoir supprimer les petites particules au centre du reservoir, qui vont se superposer avec la grosse particule ; prend en entrée la position des centres des petites particules, les vitesses de chacune des petites particules, les rayons de chacune des petites particules, les masses de chacune des petites particules, la position de la grosse particule et son rayon. 
    index_list=[] # On initialise la liste qui va contenir les index des particules que l'on doit supprimer
    points_position=config_pos.tolist() # Conversion du array des positions en liste (pour plus de simplicité dans la manipulation : pop, fonction index...)
    point_vel=config_vel.tolist() # idem avec le array des vitesses
    point_sigma=config_sigma.tolist() # idem avec le array des rayons
    point_masse= config_masse.tolist() # idem avec le array des masses
    sigma = point_sigma[0] # les rayons de toutes les petites particules étant identique, on le stock une fois pour toute dans la variable sigma
    for point in points_position[0:len(points_position)-1]: # on parcours toutes les particules, à l'exception de la dernière qui est en fait la grosse particule que l'on ne souhaite pas toucher
        if (point[0]-center_pos[0])**2+(point[1]-center_pos[1])**2<=(big_part_size+sigma)**2: # si la petite particule considérée se recouvre avec la grosse
            index = points_position.index(point) # on récupère son indice
            points_position.pop(index) # on supprime cette particule de la liste des positions
            point_vel.pop(index) # on supprime cette particule de la liste des vitesses
            point_sigma.pop(index) # on supprime cette particule de la liste des rayons
            point_masse.pop(index) # on supprime cette particule de la liste des masses
            index_list.append(index) # on stock l'index de cette particule dans la liste index
    final_pos=np.array(points_position) # on re-transforme la nouvelle liste des positions (dans laquelle on a retiré toutes les petites particules génantes) en array (car le reste des codes manipulent uniquement des arrays)
    final_vel=np.array(point_vel) # on re-transforme la nouvelle liste des vitesses (dans laquelle on a retiré toutes les petites particules génantes) en array (car le reste des codes manipulent uniquement des arrays)
    final_sigma=np.array(point_sigma) # on re-transforme la nouvelle liste des rayons (dans laquelle on a retiré toutes les petites particules génantes) en array (car le reste des codes manipulent uniquement des arrays)
    final_masse=np.array(point_masse) # on re-transforme la nouvelle liste des masses (dans laquelle on a retiré toutes les petites particules génantes) en array (car le reste des codes manipulent uniquement des arrays)
    return (final_pos,final_vel,final_sigma,final_masse,index_list) # on retourne tout les tableaux, ainsi que la liste des index supprimés

class Brownian(Simul): # On crée ici une sous classe de Simul dans laquelle on va pouvoir introduire une nouvelle grande particule
    def __init__(self, sample_time, sigma,N,L,m,M,big_part_size): # on va venir modifier la fonction __init__ de la classe simul, en y ajoutant de nouveaux paramètres en entrée : la masse de la grande particule M et sa taille big_part_size 
        super().__init__(sample_time, sigma,N,L,m) # On récupère les différentes grandeurs calculées dans la fonction __init__ de simul
        big_parti_pos = [self.taille/2,self.taille/2] # On définit la position initiale de la grande particule, au centre du réservoir
        big_parti_vel = [10**(-6),10**(-6)] # On définit la vitesse initiale de la grande particule. Cette dernière n'a pas de vitesse intrinsèque. Pour des soucis numériques (on divisera plus tard par big_parti_vel --> on veut éviter de diviser par 0), on met 10^-6 au lieu de 0 
        self.position = np.vstack([self.position,big_parti_pos]) # On rajoute à la liste des positions issue de l'init de simul une ligne pour la position de la grande particule
        self._velocity = np.vstack([self._velocity,big_parti_vel]) # On rajoute à la matrice des vitesses issue de l'init de simul une ligne pour la vitesse de la grande particule
        self.sigma = np.append(self.sigma,big_part_size)  # On rajoute à la liste des rayons issue de l'init de simul le rayon de la grande particule 
        self.masse = np.append(self.masse,M) # On rajoute à la liste des masses issue de l'init de simul la masse de la grande particule
        self.position,self._velocity,self.sigma,self.masse,index_list = remove_center_particles(self.position,self._velocity,self.sigma,self.masse,big_parti_pos,big_part_size) # on enlève les petites particules au centre du réservoir qui se recouvre avec la grande à l'aide de la fonction remove_center_particles définit au dessus de la sous-classe Brownian
        self._i, self._j = np.triu_indices(self.position.shape[0], k=1) # Renvoie une liste de deux tableaux contenant les indices des particules, et qui une fois prient ensembles donnent la totalité des coupleq de particules. Par exemple, si on a 4 particules, renvoie [array([0,0,0,1,1,2]),array([1,2,3,2,3,3])]. On le définit bien à la fin afin que la grosse particule soit compté dedans
        

        

