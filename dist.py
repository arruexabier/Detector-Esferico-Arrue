import matplotlib.pyplot as plt
import numpy as np

n=200000

archivos1=["0","A","B","C","D"]
archivos2=["0Mesh09","AMesh09","BMesh09","CMesh09","DMesh09"]
archivos3=["0Mesh08","AMesh08","BMesh08","CMesh08","DMesh08"]
casos=[archivos1,archivos2,archivos3]



for archivos in casos:
	legenda=[]
	for letra in archivos:
		hits=[]
		pos=[]
		count=0
		k=0
		z=0.5

		file=open("Caso"+letra+"_nt_Hits.csv","r")

		for i in range(8):
			next(file)
			

		for line in file:
			line_list=line.split(",")
			number=int(line_list[0])
			if number>=k:
				k=number
				count+=1
			else:
				hits.append(count/n)
				count=0
				k=0
				pos.append(z)
				z+=0.1
		print(letra,hits[0],hits[-1])
		plt.plot(pos,hits)
		legenda.append("Caso"+letra+": "+str(hits[0]))
		file.close()
	plt.xlabel("Distancia PMT-Emisión(cm)", fontsize=10)
	plt.ylabel("Probabilidad Detección", fontsize=10)
	plt.xlim([0.5,1.5])
	plt.ylim([0,1])

	if archivos==archivos1:
		def f(x):	
			r=3.2
			return (1-x/(np.sqrt(r*r+x*x)))/2.0
		y=[]
		for i in pos:
			y.append(f(i))
		#plt.plot(pos,y,"--")
		#legenda.append("Angulo sólido")

	plt.legend(legenda,loc="center left",bbox_to_anchor=(1, 0.5))
	plt.tight_layout()
	plt.yticks([0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1])
	plt.show()
	

	#ENERGY
	for letra in archivos:
		hits=[]
		pos=[]
		count=0
		k=0
		z=0.5

		file=open("Caso"+letra+"_nt_Hits.csv","r")

		for i in range(8):
			next(file)
			

		for line in file:
			line_list=line.split(",")
			number=int(line_list[0])
			if number>=k:
				k=number
				count+=1
			else:
				hits.append(3*13.7/(17*count/n))
				count=0
				k=0
				pos.append(z)
				z+=0.1
		print(letra,hits[0],hits[-1])
		plt.plot(pos,hits)
		file.close()

	plt.legend(["Caso0","CasoA","CasoB","CasoC","CasoD"],loc="center left",bbox_to_anchor=(1, 0.5))
	plt.xlabel("Distancia PMT-Emisión(cm)", fontsize=10)
	plt.ylabel("Energia (eV)", fontsize=10)
	plt.axhline(y = 13.7, color = 'y', linestyle = '--')
	plt.xlim([0.5,1.5])
	plt.tight_layout()
	plt.show()
