from src.brownian_motion import Brownian
from src.animatesimul import AnimateSimul


def main():

    simulation = Brownian(sample_time=0.015, sigma=0.15, N=36, L=5,m=1,M=10,big_part_size=0.6) # sample_time est le temps d'échantillonnage sigma,N et m sont respectivement les rayons des petites particules, leurs nombre et leurs masse. L est la taille de la boîte. M est la masse de la grande particule et big_part_size son rayon.
    
    print(simulation.__doc__)  # print the documentation from the class

    animate = AnimateSimul(simulation)
    animate.go(nframes=1000) # On choisit ici le nombre de frames.
    print(simulation)  #  print last configuration to screen

if __name__ == '__main__':
    main()
    

