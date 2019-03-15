import sqlite3

conn = sqlite3.connect('agalma-siphonophora-20180401_reduced.sqlite')
c = conn.cursor()

species16S = [ "Lilyopsis fluoracantha", "Rudjakovia sp", "Bargmannia lata", "Physonect sp", "Erenna richardi", "Resomia ornicephala", "Physophora gilmeri" ]
species18S = [ "Lychnagalma ultricularia", "Frillagalma vityazi", "Rudjakovia sp", "Bargmannia lata", "Apolemia rubriversa", "Erenna richardi", "Resomia ornicephala", "Physophora gilmeri", "Thermopalia taraxaca", "Marrus claudanielis" ]

print "\n 16S"
for species in species16S:
	export16S = open("%s_16S.fasta" % species.replace(" ", "_"), "w")
	c.execute("select c.species, m.genome_type, m.molecule_type, s.* from agalma_sequences as s join agalma_models as m on s.model_id=m.id join catalog as c on c.id=m.catalog_id where m.molecule_type='S' and m.genome_type='M' and c.species='%s';" % species)
	out16S = c.fetchall()
	print species
	print len(out16S)
	export16S.writelines([">{0}\n".format(item) for item in out16S])
	export16S.close()

c = conn.cursor()

print "\n 18S"
for species in species18S:
	export18S = open("%s_18S.fasta" % species.replace(" ", "_"), "w")
	c.execute("select c.species, m.genome_type, m.molecule_type, s.* from agalma_sequences as s join agalma_models as m on s.model_id=m.id join catalog as c on c.id=m.catalog_id where m.molecule_type='S' and m.genome_type='N' and c.species='%s';" % species)
	out18S = c.fetchall()
	print species
	print len(out18S)
	export18S.writelines([">{0}\n".format(item) for item in out18S])
	export18S.close()

conn.commit()
conn.close()