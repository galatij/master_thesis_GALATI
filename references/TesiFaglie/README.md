# Some Bibliography #

 **Libro su contact problem for linear elasticity**
Questo libro lo do solo come referenza generale. Il capitolo sul probelma di Signorini puo` servire a capire il probelma di base. Poi tratta anche Tresca e COulamb, ma forse vale la pena tenerlo sono come referenza. E` comunque uno dei libri piu` citati sull'argomento.

- `Kikuchi, Oden - Contact Probelms in Elasticity A study of variational inequalities and Finite Element Methods.pdf`

 **Book sul metodo di Nitsche** Anche qui lo do solo pre referenza, per
 avere una idea di cosa sia il metoto. Il libro diventa subito molto
 tecnico, quindi se non  interessato ad approfondire guardi solo la
 parte iniziale dove spiegano il classico metodo di Nitsche e come si
 estende pre probelmi di interfaccia per il caso del Laplaciano.
 
- `Burman and Zunino - 2011 - Numerical Approximation of Large Contrast Problems.pdf`

 **Articolo di review su problemi di contatto**
 Un altro lavoro che pu; aiutare a capire il probelma. Si concentra molto sui LAgrange MUltipliers e ignora le altre tecniche (nonostante il titolo). Di nuovo, lo terrei solo come referenza èer la èprima parte, dove il probelma è ben presentato. Parla un pò anche dei metodi di risoluzione basati su Newton semismooth, ma ci sono riferimenti migliori sull'argomento.
 
 - `Wohlmuth - Variationally consistent discretization schemes and numerical algorithms for contact problems.pdf`

 **Alcune Note**
Questa è una nota che avevo scritto per buttare giù le idee per il lavoro che poi abbiamo pubblicato. SOno delle note, qundi da prendere con le pinze, ma gliele do perché contengono anche il probelma di poroelasticità. Anche se nella tesi non andremo a quel livello (al limite aggiungeremo una pressione come data) è importante avere coscenza del contesto generale. Inoltre ci sono alcuni dettagli che per ragioni di spazio non sono contenuti nel paper.

- `NoteFaultAct.pdf`


 **Metodi basati su Lagrange multipliers (Padovani)**
Il primo è più generale, il secondo si concentra su probelma con fratture. Trattano anche il caso di frature che si aprono. Possiamo farlo anche noi, ma in una prima battuta mi limiterei a solo scorrimento.

- `JCP_LagrFaults_FraFerJanTea16.pdf`
- `Franceschini-Algebraically stabilized Lagrange multiplier method for frictional contact mechanics with hydraulically active fractures-2020-Compu.pdf`

 **Metodi basati su penalizzazioen consistente (Nitsche)**
Il primo tratta probelmi di elasticità con contatto ed è relativamente introduttivo. Ha anche un pò di analisi teorice se le interessa. Il secondo si concentra sui metodi di risoluzione. Non dobbiamo necessariamente usare i loro schemi di risoluzione, possiamo sbizzarrirci noi. E poi loro considerano molte varianti, noi ne sceglieremo una o due. Il terzo si rifereisce a un problema poromeccanico. Serve sopretutto per rivedere la formulazione attraverso una faglia discretizzata come una linea-superficie. L'ultimo è quello già fornito, con i confronti.

- `Chouly, Hild - A Nitsche-based method for unilateral contact problems Numerical analysis.pdf`
- `Chouly, Hild, Renard - A Nitsche finite element method for dynamic contact 1. Space semi-discretization and time-marching schemes.pdf`
- `Beaude et al. - Mixed and Nitsche's discretizations of Coulomb fri.pdf`
- `Chouly et al. - 2023 - Lagrangian and Nitsche methods for frictional cont.pdf
`
 **Paper Controllo**
Sopratutto per vedere la nostra notazione, diversa da quella dei Padovani.

- `Cerroni et al. - 2021 - A control problem approach to Coulomb friction.pdf`

 **Risoluzione problema non lineare**
 Qui si dà una descrizione di metodi di Newton generalizzato sia per il probelma risolto con i moltiplicatori che secondo il metodo di Nitsche.
 Presenta diverse varianti, noi vedremo di sceglierne una.
 
 - `Renard - Generalized Newton's methods for the approximation and resolution of frictional contact problems in elasticity.pdf`

