import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

## Length of the Number Piece in tenth of seconds
Nlen = 3000
## Number of iterations used for the calculations
N=10000
## Future time delay
tau=10

## List of the forte numbers for the pitch-classes. Is used for Y ticks labels
## Choose one or the other depending on the number of normal forms chosen
## in the C program.
'''
forte_names=['0-1','1-1','2-1','2-2','2-3','2-4','2-5','2-6','3-1','3-2A','3-2B',
            '3-3A','3-3B','3-4A','3-4B','3-5A','3-5B','3-6','3-7A','3-7B',
            '3-8A','3-8B','3-9','3-10','3-11A','3-11B','3-12','4-1','4-2A',
            '4-2B','4-3','4-4A','4-4B','4-5A','4-5B','4-6','4-7','4-8','4-9',
            '4-10','4-11A','4-11B','4-12A','4-12B','4-13A','4-13B','4-14A',
            '4-14B','4-z15A','4-z15B','4-16A','4-16B','4-17','4-18A','4-18B',
            '4-19A','4-19B','4-20','4-21','4-22A','4-22B','4-23','4-24','4-25',
            '4-26','4-27A','4-27B','4-28','4-z29A','4-z29B']
'''
forte_names=['0-1','1-1','2-1','2-2','2-3','2-4','2-5','2-6','3-1','3-2','3-3',
            '3-4','3-5','3-6','3-7','3-8','3-9','3-10','3-11','3-12','4-1',
            '4-2','4-3','4-4','4-5','4-6','4-7','4-8','4-9','4-10','4-11',
            '4-12','4-13','4-14','4-z15','4-16','4-17','4-18','4-19','4-20',
            '4-21','4-22','4-23','4-24','4-25','4-26','4-27','4-28','4-z29']
num_normalForms = len(forte_names)

## Data matrices
matS = np.zeros((Nlen,num_normalForms))
scaledMatS = np.zeros((Nlen,num_normalForms))
matCondS = np.zeros((Nlen-tau,num_normalForms,num_normalForms))
H = np.zeros((Nlen))
partH = np.zeros((Nlen-tau,num_normalForms))
futH = np.zeros((Nlen-tau))
mutH = np.zeros((Nlen-tau))
## the scale used to replot the probabilities
logScale = 0.25

#########################################################
### Now loading the data
#########################################################

##########
## Reading the file with instant probabilities
## matS[i,j] is P(X(i)=j)

theFile = open('outputS.txt','r')

for i in range(0,Nlen):
	for j in range(0,num_normalForms):
		matS[i,j]=float(theFile.readline())
		scaledMatS[i,j]=1.+logScale*np.log10(matS[i,j]+1./N);

scaledMatS = scaledMatS.T

#####
## Calculating the entropy H[i]=H(X(i))

for i in range(0,Nlen):
    for j in range(0,num_normalForms):
        if matS[i,j]>0:
            H[i]=H[i]-matS[i,j]*np.log2(matS[i,j])

##########
## Reading the file with conditional probabilities
## matCondS[i,j,k] is P(X(i+tau)=j | X(i)=k)

theFile = open('outputCondS.txt','r')

for i in range(0,Nlen-tau):
	for j in range(0,num_normalForms):
         for k in range(0,num_normalForms):
             matCondS[i,j,k]=float(theFile.readline())

#####
## Calculating the partial entropy partH[i,j]=H(X(i+tau)|X(i)=j)

for i in range(0,Nlen-tau):
	for j in range(0,num_normalForms):
         for k in range(0,num_normalForms):
             if matCondS[i,k,j]>0:
                 partH[i,j] = partH[i,j]-matCondS[i,k,j]*np.log2(matCondS[i,k,j])

#####
## Calculating the future entropy futH[i]=H(X(i+tau)|X(i))

for i in range(0,Nlen-tau):
        for j in range(0,num_normalForms):
                futH[i] = futH[i]+matS[i,j]*partH[i,j]

mutH = H[tau:Nlen]-futH
partH=partH.T

#########################################################
### Now displaying the results
#########################################################

### Definition of the appropriate colormap


### This is black-violet-red-yellow-white
cdict = {'red': ((0.0, 0.0, 0.0),
				 (0.25, 0.328, 0.328),
				 (0.5, 0.945, 0.945),
				 (0.75, 0.855, 0.855),
				 (1.0, 1.0, 1.0)),
		 'green': ((0.0, 0.0, 0.0),
				 (0.25, 0.0, 0.0),
				 (0.5, 0.105, 0.105),
				 (0.75, 0.792, 0.792),
				 (1.0, 1.0, 1.0)),
		 'blue': ((0.0, 0.0, 0.0),
				 (0.25, 0.629, 0.629),
				 (0.5, 0.0, 0.0),
				 (0.75, 0.0, 0.0),
				 (1.0, 1.0, 1.0))}

'''
### This is black-blue-green-yellow-white
cdict = {'red': ((0.0, 0.0, 0.0),
				 (0.25, 8./255., 8./255.),
				 (0.5, 0., 0.),
				 (0.75, 199./255., 199./255.),
				 (1.0, 1.0, 1.0)),
		 'green': ((0.0, 0.0, 0.0),
				 (0.25, 25./255., 25./255.),
				 (0.5, 138./255., 138./255.),
				 (0.75, 187./255., 187./255.),
				 (1.0, 1.0, 1.0)),
		 'blue': ((0.0, 0.0, 0.0),
				 (0.25, 189./255., 189./255.),
				 (0.5, 36./255., 36./255.),
				 (0.75, 0.0, 0.0),
				 (1.0, 1.0, 1.0))}

### This is darkpurple-blue-turquoise-white
cdict = {'red': ((0.0, 0.0, 0.0),
				 (0.25, 72./255., 72./255.),
				 (0.5, 60./255., 60./255.),
				 (0.75, 111./255., 111./255.),
				 (1.0, 1.0, 1.0)),
		 'green': ((0.0, 0.0, 0.0),
				 (0.25, 0./255., 0./255.),
				 (0.5, 112./255., 112./255.),
				 (0.75, 200./255., 200./255.),
				 (1.0, 1.0, 1.0)),
		 'blue': ((0.0, 0.0, 0.0),
				 (0.25, 173./255., 173./255.),
				 (0.5, 229./255., 229./255.),
				 (0.75, 204./255., 204./255.),
				 (1.0, 1.0, 1.0))}
'''
my_cmap = LinearSegmentedColormap('my_colormap',cdict,256)

######
## Image picture of matS

theFig = plt.figure(figsize=(14,9))
theImage = plt.imshow(scaledMatS,aspect='auto',interpolation='none',cmap=my_cmap)
theImage.set_clim(0.0,1.0)
plt.xlabel('Time (tenths of seconds)',fontsize=10)
plt.xticks(fontsize=10)
plt.yticks(range(0,num_normalForms),forte_names,fontsize=10)

### Everything we need for the colorbar

colorbarTicks=np.arange(0.0,1.1,0.1)
colorbarTicksLabels = []
for i in colorbarTicks:
    colorbarTicksLabels.append("%.1e" % 10.0**((i-1.0)/logScale))
##colorbarTicksLabels[0]="%.2e" % 0.0

theCB = plt.colorbar()
theCB.set_ticks(colorbarTicks)
theCB.set_ticklabels(colorbarTicksLabels)
theCB.ax.tick_params(labelsize=10)

######
## Now we display the conditional entropy image

theFig2 = plt.figure(figsize=(10,10))
theImage2 = plt.imshow(partH,aspect='auto',interpolation='none',cmap=my_cmap)
plt.xlabel('Time (tenths of seconds)',fontsize=10)
plt.xticks(fontsize=10)
plt.yticks(range(0,num_normalForms),forte_names,fontsize=10)
theCB2 = plt.colorbar()
theCB2.ax.tick_params(labelsize=10)

######
## Now we display the entropy plot
## Black: Entropy H(X(i))
## Gray: future entropy H(X(i+tau)|X(i))
## Dashed: mutual entropy I(X(i+tau),X(i))

theFig3 = plt.figure()
thePlot = plt.plot(range(0,3000),H,'-k')
plt.ylim([0,6])
plt.plot(range(0,3000-tau),futH,color='0.5')
plt.plot(range(0,3000-tau),mutH,'--k')
plt.xlabel('Time (tenths of seconds)',fontsize=16)
plt.ylabel('Entropy (in bits)',fontsize=16)
plt.grid()

## Done
plt.show()
