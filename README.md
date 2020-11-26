# projet_Maths

Codage en C.

On utilise stb_image pour charger une image et la convertir en niveaux de gris.

- Transformée de fourrier directe discrète sur une image
  
  double boucle
- Transformée de fourrier rapide discrète 1D

  ATENTION: Le signal doit avoir 2^n échantillons sinon la décomposition/recomposition ne peut fonctionner.

  un tableau pour les valeurs & un tableau pour les index

  On mélange l'entrée en inversant les bits: 110 -> 011 et on stocke les index

  On va avoir une fonction pour faire une "mini TF" entre les recompositions (prend un tableau, deux indices et retourne un double)

  On va utiliser la méthode butterfly (a voir si il faut la faire récursive)

  On stocke les résultats dans un tableau de sortie (double)
  
- Transformée de fourrier rapide discrète 2D

  On peut la séparer en deux FFT 1D, c'est a dire FFT(FFT(image))

- Transformée de fourrier inverse rapide 1D

- transformée de fourrier inverse rapide 2D
