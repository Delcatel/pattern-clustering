# pattern-clustering

Ceci est un projet de recherche de 1ère année à l'Ecole Nationale des Ponts et Chaussées. Il est attaché au cours d'ouverture de traitement du signal et pose plus particulièrement le problème de l'estimation aveugle de canal.

Les données d'entrée sont obtenu par la génération aléatoire de points dans un domaine borné du plan. Puis par l'implémentation d'algorithme inspiré par le k-means (ou k-product), on "clusterise" en imposant une condition géométrique sur les modes (ou centroïdes) des différents clusters :
  * appartenir à un cercle
  * former un triangle
  * former un carré
