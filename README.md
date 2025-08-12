# HeatMaster

1. Działanie programu

Program umożliwia numeryczną analize nieustalonego pola temperatury.
W tym celu program wykorzystuje równania cieplne modeli otrzymane dzięki metodzie różnic skończonych.
Metoda obliczeniowa jest zależna od wprowadzonych danych wejściowych. Jest to metoda jawna i iteracyjna.
Z tego też powodu wprowadzony układ może być niestabilny, wtedy można próbować rozwiązać ten problem zmniejszająć krok czasowy analizy.
Należy pamiętać, że metoda jawna jest 2 razy bardziej stabilna, niż metoda iteracyjna tzn. że przy dwukrotnie większym kroku czasowym obiekt powinien być stabilny, lecz jest to kosztem dokładności obliczeń.

W programie są dostępne tylko dwa modele:
- ogrzewanie podłogowe - z założeniem, że kanaliki grzewcze są równomiernie rozłożone i są w środku płyty betonowej. Analizowany jest tylko fragment podłogi, który może być ekstrapolowany na całą resztę.
- płaski szynoprzewód - Analiza temperaturowa w tym modelu ogranicza się do prądu stałego, ponieważ nie uwzględnia zjawiska naskórkowości. Analizowany fragment jest ograniczony tylko dla górnej-prawej ćwiartki, która może być ekstrapolowane na cała resztę szynoprzewodu.

Oprawa graficzna została wykonana przy użyciu biblioteki wxWidgets.
Program dokonuje obliczeń równolegle na karcie graficznej NVIDIA przy użyciu biblioteki CUDA.

2. Obsługa programu

Po uruchomieniu się okna programu z lewej górnej strony znajduje wybór modelu. Po prawej stronie stronie tego obiektu znajduje się okno pozwalające wybrać, czy obliczenia powinny być wykonywane przy pomocy macierzy gęstych, czy rzadkich pozwalających zaoszczędzić pamięć operacyjną. Na tej samej wysokości na prawym brzegu znajduje się przycisk "Uruchom" rozpoczynający przeprowadzanie obliczeń.
Aby wprowadzić parametry wybranego modelu i analizy należy na górnym pasku kliknąć "Parametry>Wprowadź dane".
Po przeprowadzeniu analizy na dolnym ekranie pojawiają się kolory od barwy niebieskiej do czerwonej. Symbolizują one temperaturę w skali 10 do 100 stopni Celsjusza.
Nie dodano do programu dostosowywującej się do sytuacji skali temperatury, ani legendy ułatwiającej odczyt temperatury obiektu. Choć są ważne to jednak zostały pominięte, ponieważ nie były obiektem badań.
Aby zmienić wyświetlaną temperaturę w danej chwili czasowej wystarczy przesuwać suwakiem znajdującym się nad wyświetleniem kolorów.
