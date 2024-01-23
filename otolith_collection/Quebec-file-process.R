##
## "//dcqcimlna01a/BD_Peches","/Releves_Poissons_de_Fond_et_Crevette/Extractions_Oracle/RelevésRecherche/Formats d'extraction"
##
process.que.file.fct <- function(in.str, file.type){
    switch(file.type,
           sets={
               df <- data.frame(
                   nbpc=trimws(substring(in.str, 1,2), which="left"),
                   no.rel=trimws(substring(in.str, 3,5), which="left"),
                   trait=trimws(substring(in.str, 6,8), which="left"),
                   date6=substring(in.str, 9,14),
                   ttow=substring(in.str, 15,15),
                   no.str=substring(in.str, 16,18),
                   n.unit=substring(in.str, 19,21),
                   nafo=substring(in.str, 22,23),
                   light=substring(in.str, 24,24),
                   nebulo=substring(in.str, 25,25),
                   dir.ven=substring(in.str, 26,28),
                   for.ven=substring(in.str, 29,29),
                   weather=substring(in.str, 30,31),
                   mer=substring(in.str, 32,32),
                   pres.atm=substring(in.str, 33,35),
                   air.t.sec=substring(in.str, 36,39),
                   air.t.hum=substring(in.str, 40,43),
                   per.vag=substring(in.str, 44,45),
                   hau.vag=substring(in.str, 46,47),
                   dir.hou=substring(in.str, 48,50),
                   per.hou=substring(in.str, 51,52),
                   hau.hou=substring(in.str, 53,54),
                   t.hr=substring(in.str, 55,55),
                   fus.hor=substring(in.str, 56,56),
                   hr.debut=as.numeric(substring(in.str, 57,58)),
                   min.debut=substring(in.str, 59,60),
                   hr.fin=as.numeric(substring(in.str, 61,62)),
                   min.fin=substring(in.str, 63,64),
                   vit=substring(in.str, 65,67),
                   dir.tow=substring(in.str, 68,70),
                   dist.towed=as.numeric(trimws(substring(in.str, 71,73), which="left")),
                   funes=substring(in.str, 74,76),
                   resul=trimws(substring(in.str, 77,78), which="left"),
                   prof.min=substring(in.str, 79,82),
                   prof.max=substring(in.str, 83,86),
                   prof.moy=substring(in.str, 87,90),
                   prof.bot=substring(in.str, 91,94),
                   prof.fis=substring(in.str, 95,98),
                   temp.sur=substring(in.str, 99,101),
                   temp.fond=substring(in.str, 102,104),
                   lat.debut=substring(in.str, 105,109),
                   lon.debut=substring(in.str, 110,114),
                   lat.fin=substring(in.str, 115,119),
                   lon.fin=substring(in.str, 120,124),
                   mt.pos=substring(in.str, 125,125),
                   engin=substring(in.str, 126,127),
                   sal.fond=substring(in.str, 128,132),
                   oxygene=substring(in.str, 133,139), stringsAsFactors=F
               )

			   df$lon.debut.dec <- -1* (
                   as.numeric(substring(df$lon.debut,1,2)) +
                   ((as.numeric(substring(df$lon.debut,3,4))+(0.1*as.numeric(substring(df$lon.debut,5,5))))/60)
               )

               df$lon.fin.dec <- -1 * (
                   as.numeric(substring(df$lon.fin,1,2)) +
                   ((as.numeric(substring(df$lon.fin,3,4))+(0.1*as.numeric(substring(df$lon.fin,5,5))))/60)
               )

               df$lat.debut.dec <- (
                   as.numeric(substring(df$lat.debut,1,2)) +
                   ((as.numeric(substring(df$lat.debut,3,4))+(0.1*as.numeric(substring(df$lat.debut,5,5))))/60)
               )

               df$lat.fin.dec <- (
                   as.numeric(substring(df$lat.fin,1,2)) +
                   ((as.numeric(substring(df$lat.fin,3,4))+(0.1*as.numeric(substring(df$lat.fin,5,5))))/60)
               )

           },
           catch={

               df <- data.frame(
                   nbpc=trimws(substring(in.str, 1,2), which="left"),
                   no.rel=trimws(substring(in.str, 3,5), which="left"),
                   trait=trimws(substring(in.str, 6,8), which="left"),
                   date6=as.numeric(substring(in.str, 9,14)),
                   Esp=as.numeric(substring(in.str, 18,21)),
                   NbIndiv=as.numeric(substring(in.str, 22,27)),
                   Pds.Capt=as.numeric(substring(in.str, 28,39))/1000000,
                   Pds.Ech=as.numeric(substring(in.str, 40,51))/1000000
                   )
           },
		   lengthSTRAP={
               df <- data.frame(
                   nbpc=trimws(substring(in.str, 1,2), which="left"),
                   no.rel=trimws(substring(in.str, 3,5), which="left"),
                   no.str=trimws(substring(in.str, 6,8), which="left"),
                   trait=trimws(substring(in.str, 9,11), which="left"),
                   date6=as.numeric(substring(in.str, 12,17)),
                   division=substring(in.str, 18,19),
                   ttow=substring(in.str, 20,20),
                   Esp=as.numeric(substring(in.str, 21,24)),
                   Nb.Mes=substring(in.str, 25,28),
                   ratio=substring(in.str, 29,30),
                   sex=substring(in.str, 31,31),
                   group=substring(in.str, 35,35),
                   lf=substring(in.str, 39,356), stringsAsFactors=F
                   )
           },
		   length={
			df <- data.frame(
				nbpc=trimws(substring(in.str, 1,2), which="left"),
				no.rel=trimws(substring(in.str, 3,5), which="left"),
				no.station=trimws(substring(in.str, 6,8), which="left"),
				date6=as.numeric(substring(in.str, 9,14)),
				Esp=as.numeric(substring(in.str, 15,18)),
				sex=substring(in.str, 19,19),
				length.min=substring(in.str, 20,22),
				length.max=substring(in.str, 23,25),
				category=substring(in.str, 26,28),
				Nb.Mes=substring(in.str, 29,31)
				)
           },
			bio={
               df <- data.frame(
                   nbpc=trimws(substring(in.str, 1,2), which="left"),
                   no.rel=trimws(substring(in.str, 3,5), which="left"),
                   trait=trimws(substring(in.str, 6,8), which="left"),
                   date6=as.numeric(substring(in.str, 9,14)),
                   Esp=as.numeric(substring(in.str, 15,18)),
                   Specimen=as.numeric(substring(in.str, 19,23)),
                   Longueur=as.numeric(substring(in.str, 24, 27)),
                   sex=substring(in.str,41,41),
                   Poids=as.numeric(substring(in.str, 49, 56)),
                   Longueur2=as.numeric(substring(in.str, 253, 256))
                   )
           },
		   prelevement={
		   df <- data.frame(
                   nbpc=trimws(substring(in.str, 1,2), which="left"),
                   no.rel=trimws(substring(in.str, 3,5), which="left"),
                   trait=trimws(substring(in.str, 6,8), which="left"),
                   date6=as.numeric(substring(in.str, 9,14)),
                   Esp=as.numeric(substring(in.str, 15,18)),
                   Specimen=as.numeric(substring(in.str, 19,23)),
                   Sexe=as.numeric(substring(in.str, 24,24)),
				   Type.prelevement=as.numeric(substring(in.str, 26,27)),
				   Type.demande=substring(in.str, 29,30),
				   Tri=substring(in.str, 31,32),
				   Categorie=substring(in.str, 33,34)
                   )
		   },
           stop("File type must be sets, catch, lengthSTRAP, length, bio or prelevement")
           )
return(df)
}

