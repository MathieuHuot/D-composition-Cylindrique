# Décomposition-Cylindrique
L'ensemble des algorithmes formant la décomposition cylindrique dans le but de l'application aux PolITA

Authors : Mathieu Huot
          Rémy Garnier

Context : Dans le cadre de notre stage de L3 pour la licence de mathématique à l'ENS Cachan, nous étudions au                    laboratoire du LSV l'implémentation d'un algorithme résolvant la logique du premier ordre des réels dans le            cas de réels algébriques.

Voici dans une première partie l'avancement de notre stage .
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


5) juin 1-15
~~~~~~~~~~~~

-système de plot de résultats du lifting
-mise en contact avec un ingénieur INRIA pour l'interfaçage avec Cosy-Verif
-introduction à Lupa
-ajout d'une simplification randomizée dans elim
-mise en contact avec l'ingénieur du LSV pour la necessité de puissance de calcul
-parallélization de elim et line partition
-stockage de root-coding pour limiter redondances de calculs


6) juin 16-30
~~~~~~~~~~~~~

-dérécursivation de lifting et gain en mémoire

II) Difficultés surmontées
==========================
