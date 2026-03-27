------------------------------------
1. Quellcode
------------------------------------
-Das Projekt wurde in C++ geschrieben
-die .cpp Datein sind in "src" zu finden
-die zugehörigen .h Dataein in "include"
-im main.cpp befindet sich die main methode
  in der die verschiedenen Simultionen zu den jeweiligen Aufgaben vorzufinden sind
  -durch commenting lassen sich die Aufgaben einzeln ausführen
  -wie man die programmierten Integratoren benutzt steht dort auch als Kommentar dabei 
  -per default ist jede Aufgabe zunächst auskommentiert außer die Initialisierung der Integratoren.
-bisjetzt haben wir den konstanten und dynamischen Zeitschritt umgesetzt. Wir haben jedoch sowieso nur konstante Zeitschritte benutzt. 
------------------------------------
2. Wichtiges
------------------------------------
Bei der Generierung unser Files hat sich ein Bug eingeschlichen:
-In der Datei im Header steht "... E	p_x	p_y	p_x	p_z  ..."
-Dies führt zu 18 columns, jedoch haben gibt es nur 17 die tastächliche Werte beinhalten.
-Wenn man in python zum Einlesen der Dateiten Pandas (pd.read_table(...)) führt die zu einem column a am Ende der Leer ist.
-ab dem zweiten p_x sind alle columns um eins nach links zu verschieben. (z.B Die Daten von p_z sind eigentlich die von l_x, die Daten von |j| sind die von |e|)
-Dies ist uns leider erst später nach den Simulatinen aufgefallen...
-----------------------------------------------
3. Kompilierung und Starten des Programms
-----------------------------------------------
