# Décomposition-Cylindrique
L'ensemble des algorithmes formant la décomposition cylindrique dans le but de l'application aux PolITA

Authors : Rémy Garnier
          Mathieu Huot
          
Context : Dans le cadre de notre stage de L3 pour la licence de mathématique à l'ENS Cachan, nous étudions au                    laboratoire du LSV l'implémentation d'un algorithme résolvant la logique du premier ordre des réels dans le            cas de réels algébriques.

Voici dans une première partie l'avancement de notre stage.
Dans une seconde partie les difficultés que nous avons eues avec l'implémentation.

I) Avancement du stage
======================

1) janvier février
~~~~~~~~~~~~~~~~~~

-auto formation à Python
-auto formation à Sage
-cours et réunion avec encadrants pour comprendre l'article
-compréhension des preuves et des algorithmes
-correction de coquilles et imprécisions dans l'article

2) mars 
~~~~~~~

-utilisation de Sagecloud et Overleaf
-implémentation de la phase de Sign dans Q[X] (les 7 algorithmes récursifs)
-correction de coquilles et imprécisions dans l'article

3) avril
~~~~~~~~

-calculs de complexité de certains algorithmes
-implémentation de la décomposition cylindrique dans Q[X1,X2] (Elim, Lifting, LinePartition et Completing)
-utilisation de beamer latek, géogebra, wolfram alpha
-préparation de la soutenance de mi-stage

4) mai
~~~~~~

-décomposition cylindrique dans Q[X1,...,Xn]
-auto formation à Git, Github
-ajout de conversion d'anneaux, recodage de quotient
-simplification de Elim avec les pgcd
-simplification dans Line Partition et Completing avec les précalculs
-introduction aux PolITA : interets, antécédants, cadre de l'étude
-refonte du code de base : plus de commentaires et de clarté, suppression de fonctions devenues inutiles

5) juin 1-14
~~~~~~~~~~~~

-système de plot de résultats du lifting
-mise en contact avec un ingénieur INRIA pour l'interfaçage avec Cosy-Verif
-introduction à Lupa
-ajout d'une simplification randomizée dans elim
-mise en contact avec l'ingénieur du LSV pour la necessité de puissance de calcul
-parallélization de elim, line partition
-simplification dans Elim avec les pgcd
-simplification de calcul de Sign dans Lifting
-implémentation d'un dictionnaire pour l'algorithme "on the fly"
-implémentation de la classe PolITA
-amélioration générale de la qualité du code
-amélioration de SignRealization
- interface LUA : nous avons installé lupa, et essayé de faire marcher les tables. Il y avait un problème d'indicage avec les entiers Sage qui sont différents des entiers python. On a posé la question sur le forum Sage pour savoir comment contourner le problème.
- algorithme "on the fly" et algorithme d'accessibilité
- parallelisation de completing, lifting
- avancées dans les PolITA

6) juin 15-29
~~~~~~~~~~~~~


II) Difficultés surmontées
==========================

-isomorphisme en mathématiques des anneaux Q[X1][X2]...[Xn], Q[X1,...,Xn], Q[X1,...,Xn-1][Xn] mais bien entendu pas en programmation !
-découverte assez tardive de l'ordre monomial sur les variables et des bases de Gröbner pour la compréhesnion de certaines opérations d'anneaux
-la division d'un polynome par 1 fait une fraction rationnelle 
-Mais on peut forcer l'anneau dans lequel le polynome est considéré, on peut lui ajouter artificiellement +0*Xi pour le considérer dans un certain anneau
-peu de documentation et des gens soit trop peu connaisseurs soit trop connaisseurs pour nous aider
-la facilité d'utilisation de Python est contrebalancée par un manque de typage qui demande parfois un debugging assez long
-sage est parfois trop typé
 
