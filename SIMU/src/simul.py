import numpy as np

class Simul:
    """ 
    This is the prototype of the simulation code
    It moves the particles with at _velocity, using a vector notation: numpy should be used.
    """
    def __init__(self, sample_time, sigma,N,L,m): # prend en entrée l'instance de la classe, le temps d'échantillonage, le rayon des petites particules (un scalaire), le nombre de petites particules, le la longueur du coté du réservoir (carré) et la masse des petites particules (un scalaire)
        print("Simul init")
        # On initialise le système à un cristal, où les petites particules sont régulièrement espacées.
        n = int(np.ceil(np.sqrt(N))) # On définit le nombre de particules qu'il y aura sur une ligne/colonne (exemple : si N = 36, alors 6 particules par ligne et par colonnes) 
        self.taille = L # On fait de L une variable self pour pouvoir la réutiliser dans animatesimul par la suite
        absi = np.linspace(self.taille/(n+1),self.taille*n/(n+1),n)  # Séquençage en abscisses : on subdivise L en n points
        ordo = np.linspace(self.taille/(n+1),self.taille*n/(n+1),n)  # Séquençage en ordonnées : on subdivise L en n points
        position = [] # on initialise la liste des positions des petites particules
        for i in range(len(absi)):
                       for j in range(len(ordo)):
                               position.append([absi[i],ordo[j]])  # Création de la liste des positions en cristal des petites particules à t=0 : on place le centre des petites particules sur chacun des n carré points préalablement définis
        self.position=np.array(position[:N]) # Transformation de position en array et en variable globale de classe
        self._velocity = np.random.normal(size=self.position.shape) # Création de la liste des vitesses sous forme d'une variable globale de classe pour chaque petite particule : vitesse aléatoire
        a,b = self.position.shape # Taille du tableau des positions
        self.sigma = np.full(a, sigma) # Définition d'une liste de rayon, de même taille que la liste des positions que l'on remplit avec le scalaire sigma. Toutes les petites particules ont le même rayon. On peut évidemment jouer sur ce paramètre comme on veut (d'où la force d'utiliser une liste pour sigma plutôt qu'un simple scalaire).
        self.masse = np.full(a, m) # Définition d'une liste de masse, de même taille que la liste des positions que l'on remplit avec le scalaire m. Toutes les petites particules ont la même masse. Comme pour les rayons, on pourrait jouer sur cela pour avoir des petites particules de masses toutes différentes.
        self._i, self._j = np.triu_indices(self.position.shape[0], k=1) # Renvoie une liste de deux tableaux contenant les indices des particules, et qui une fois prient ensembles donnent la totalité des coupleq de particules. Par exemple, si on a 4 particules, renvoie [array([0,0,0,1,1,2]),array([1,2,3,2,3,3])]
        self._sample_time=sample_time # Temps d'échantillonage sous forme d'une variable globale de classe
        self.position_big_part = np.array([[self.taille/2,self.taille/2]]) # On place la position initiale de la grosse particule dans une liste, à laquelle on ajoutera par la suite sa position au cours du temps --> variable self car on la réutilisera dans d'autre fonction par la suite, en faisant simplement appelle à self ; n'est pas utilisée quand simulation sans grosse particule
             
    def _wall_time(self): # Doit renvoyer le temps au bout duquel il y a la première collision entre un mur et une particule, la particule concernée et le mur collisionné ; prend en entrée la variable globale de classe (instance de classe)
        time = [] # initialisation de la liste qui contiendra par la suite les temps de collision de chaque particule avec un mur
        mur=[] # initialisation de la liste qui contiendra l'information sur les murs percutés par chacune des particules 
        for k in range(self.position.shape[0]): # on considère chacune des petites particules 
            temps_first=[(self.taille-self.sigma[k]-self.position[k][0])/self._velocity[k][0],(self.taille-self.sigma[k]-self.position[k][1])/self._velocity[k][1],(self.sigma[k]-self.position[k][0])/self._velocity[k][0],(self.sigma[k]-self.position[k][1])/self._velocity[k][1]] # On calcule pour chaque particule, en tenant compte de sa position, son rayon et de sa vitesse les temps nécessaires pour qu'elle tape chacun des murs. temps_first est donc une liste de 4 temps correspondant au temps au bout desquels la kième particule tape le mur du haut, de bas, de droite et de gauche (pas nécessairement dans cet ordre). On l'utilisera par la suite pour savoir dans quel mur la particule tape en premier.
            temps=[t for t in temps_first if t>0] # De la liste temps_first, on enlève les temps négatifs qui ne sont pas physiques (ils reviendraient à avoir une vitesse négative).
            min_time=min(temps) # On prend le plus petit de ces temps, correspondant au temps au bout duquel la particule collisionnera une première fois un des murs (les autres n'ont pas d'importance, la particule changera de trajectoire après sa collision avec un premier mur
            mur_dot=temps_first.index(min_time) # On récupère l'indice auquel le minimum est atteint dans la liste temps_first. Vu la méthode de calcul de cette liste, on pourra remonter à l'information "mur tapé par la particule" (si 0 ou 2, mur vertical, si 1 ou 3, mur horizontal)
            # Conversion de l'index en 0 si le mur tapé est vertical ou 1 si le mur tapé est horizontal
            if mur_dot == 0 or mur_dot ==2 : 
                mur_dot=0 # Mur vertical
            if mur_dot == 1 or mur_dot ==3 :
                mur_dot=1 # Mur horizontal
            time.append(min_time) # Dans la liste time, on récupère les temps de collisions de chaque particule de notre système avec un mur
            mur.append(mur_dot) # Dans la liste time, on récupère l'information sur le mur tapé par chacune des particules
            min_min_time=min(time) # C'est le temps qui nous intéresse et qu'on retournera : c'est le plus petit temps au bout duquel l'une des particules touche l'un des murs. Ce temps nous permettra de tout actualiser dans la fonction md_step
            min_index=np.argmin(time) # On identifie laquelle des particules touche un mur au premier
            wall_index=mur[min_index]  # On identifie quel mur est percuté par cette particule
        return min_min_time,min_index,wall_index 

    def _pair_time(self):  # Le but de cette fonction est de retourner le plus petit temps au bout duquel il y a collision entre deux particules du système ainsi que les deux particules concernées ; prend en entrée la variable globale de classe
        r_ij = self.position[self._i] - self.position[self._j] # Vecteur (liste de liste à deux élèments) des différences des positions entre les particules ; par exemple si on prend 4 particules, r_ij ressemblerait à [[Delta_x(12),Delta_y(12)],[Delta_x(13),Delta_y(13)],[Delta_x(14),Delta_y(14)],[Delta_x(23),Delta_y(23)],[Delta_x(24),Delta_y(24)],[Delta_x(34),Delta_y(34)]] en notant (ij) les particules concernées i et j
        v_ij = self._velocity[self._i] - self._velocity[self._j] # Vecteur (liste de liste à deux élèments) des différences des vitesse des particules. Même principe que r_ij.
        # Pour chaque couple de particules, on calcule le temps qu'il faut pour qu'ils rentrent en collision, en fonction de leurs positions, vitesses et rayons. Quantitativement, ça revient à résoudre une équation du second degré et donc le calcul d'un discriminant.
        sum_sigma = self.sigma[self._i] + self.sigma[self._j] # Définition du terme de droite de l'équation quadratique
        Det = [4*(np.dot(r_ij[k],v_ij[k]))**2-4*(np.linalg.norm(v_ij[k]))**2*(np.linalg.norm(r_ij[k])**2-sum_sigma[k]**2) for k in range(len(r_ij))] # Calcul du déterminant de l'équation quadratique
        Sol = [] # Initialisation de la liste des solutions de l'équation du second degré
        for k in range(len(Det)): # Pour chaque couple de particules
            if 2*(np.dot(r_ij[k],v_ij[k])) < 0 : # on vérifie que le coefficient "b" soit négatif pour la résolution de l'équation du second ordre (comme ça -b est positif et on a pas de problème avec des solutions négatives)
                if Det[k]>=0: # Si le determinant est positif ou nul (car on le passe à la racine dans l'expressionn des solutions)
                    Solu_1 = (-2*np.dot(r_ij[k],v_ij[k]) - np.sqrt(Det[k]))/(2*np.linalg.norm(v_ij[k])**2)
                    Solu_2 = (-2*np.dot(r_ij[k],v_ij[k]) + np.sqrt(Det[k]))/(2*np.linalg.norm(v_ij[k])**2)
                    # On choisit le minimum des temps positifs 
                    if Solu_1 > 0 and Solu_2 > 0: # Si les deux solutions sont positives
                        Sol.append(min(Solu_1,Solu_2)) # On choisi la plus petites des deux (on veut le premier temps de collision)
                    elif Solu_1*Solu_2<0 : # Si les deux solutions sont de signe opposé (une positive et une négative), on prend la plus grande des deux, ie la solution positive
                        Sol.append(max(Solu_1,Solu_2)) 
                else: # Si les deux solutions sont négatives
                    Sol.append(10**5) # On place dans sol un temps très important afin d'être sur dans md_step qu'il ne soit jamais considéré comme temps de collision
            else : # si le coefficient b est négatif (pas de solution positive)
                Sol.append(10**5) # On place dans sol un temps très important pour dire qu'il n'y aurait pas de collision 
        min_time_collision,collision_index = min(Sol), np.argmin(Sol) # On prend le minimum des temps de collisions entre chaque paire de particules, et on récupère l'index dans sol de ce minimum
        # On récupère les indices des deux particules qui se colissionnent en min_time_collision
        i = self._i[collision_index] 
        j = self._j[collision_index]      
        return(min_time_collision,i,j)
        
    def md_step(self): # On souhaite pouvoir actualiser la position de chacune des particules après chaque évènement de collision (entre une particule et un mur ou entre deux particules) ; prend en entrée la variable globale de classe
        print('Simul::md_step')
        pressure = 0 # On initialise la pression à 0
        temps_residuel = self._sample_time # Au début, le temps résiduel, ie le temps pendant lequel on continue d'actualiser la position de chacune des particules, est le temps de notre sampling (temps d'échantillonage)
        temps_mur, part_mur, mur = self._wall_time() # On applique une première fois la fonction wall_time 
        temps_collision,i,j=self._pair_time() # On applique une première fois la fonction pair_time 
        while min(temps_mur,temps_collision) < temps_residuel : # Tant que le minimum entre le temps avant une collision particule/mur et le temps avant une collision particule/particule est plus petite que le temps résiduel (ie il peu encore y avoir une collision d'un des deux types pendant le temps restant) : 
            #Si le temps avant une collision particule/mur est plus petit que le temps avant une collision particule/particule, on réinitialise le système avec la donnée du temps de collision avec le mur. On fera l'inverse quand si temps collision part/part < temps collision part/mur.
            if temps_mur < temps_collision : 
                temps_residuel = temps_residuel - temps_mur # Le nouveau temps résiduel est désormais défini comme le temps restant après la collision avec le mur
                self.position = self.position + temps_mur * self._velocity # On actualise la position des particules (d=v*t)
                self._velocity[part_mur,mur] = -self._velocity[part_mur,mur] # Le sens de la vitesse s'inverse, avec conservation de la norme de la vitesse.
                pressure+=(1/self._sample_time)*self._velocity[part_mur,mur] # La pression de notre système augmente selon la relation décrite dans le fichier Pressure.ipynb
            else: # Si au contraire le temps de collision part/part est plus faible que le temps de collision part/mur
                part = [i,j] # Les deux particules entre lesquelles il y aura collision
                temps_residuel = temps_residuel - temps_collision # Le nouveau temps résiduel est désormais défini comme le temps restant après la collision entre les deux particules
                self.position = self.position + temps_collision * self._velocity # On actualise la position (d=v*t)
                vect_non_unitaire = self.position[part[1]] - self.position[part[0]] # Calcul du vecteur qui relie les deux particules qui rentrent en collision
                vect_unitaire= vect_non_unitaire/np.linalg.norm(vect_non_unitaire) # Normalisation du vecteur reliant les deux particule
                self.diff_vit = self._velocity[part[0]] - self._velocity[part[1]] # Vecteur différence de vitesse des deux particules entrant en collision
                # On actualise les vitesses de chacune des particules après leurs collision. 
                self._velocity[part[0]] = self._velocity[part[0]] - 2*(self.masse[j]/(self.masse[i]+self.masse[j]))*np.dot(self.diff_vit,vect_unitaire)*vect_unitaire # Première particule
                self._velocity[part[1]] = self._velocity[part[1]] + 2*(self.masse[i]/(self.masse[i] + self.masse[j]))* np.dot(self.diff_vit,vect_unitaire)*vect_unitaire # Seconde particule
            temps_collision,i,j = self._pair_time() # On ré-execute une fonction pair_time, sous ces nouvelles conditions, toujours dans la boucle
            temps_mur, part_mur, mur = self._wall_time() # On ré-execute une fonction wall_time, sous ces nouvelles conditions, toujours dans la boucle
            # Et on repasse dans la boucle en comparant le nouveau temps résiduel avec chacun des temps t_col_part/mur et t_col_part/part
        # En sortant de la boucle, cela signifie que le temps_residuel n'est plus suffisament grand pour permettre pour une nouvelle collision (part/mur ou part/part), on va donc simplement actualiser le système en utilisant le temps résisuel (v=dt)
        self.position = self.position + temps_residuel * self._velocity 
        # On retourne la position et la vitesse de la grande particule pour analyser son mouvement dans le bain des petites particules.
        self.position_big_part = self.position[-1] # Renvoie la position de la grosse particule en sortie de la fonction md_step (donc tous les temps d'échantillonage)
        self.velocity_big_part = self._velocity[-1] # Renvoie la vitesse de la grosse particule en sortie de la fonction md_step (donc tous les temps d'échantillonage)
        return (self.position_big_part,self.velocity_big_part)
    
    def __str__(self):   # this is used to print the position and velocity of the particles
        p = np.array2string(self.position)
        v = np.array2string(self._velocity)
        #Sol = np.array2string(self._pair_time())
        return 'pos= '+p+'\n'+'vel= '+v+'\n'+'times'
    #+'times='+sol

