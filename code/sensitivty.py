import pandas as pd
import numpy as np
from sklearn.ensemble import RandomForestRegressor
from sklearn.datasets import make_regression
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.cross_validation import cross_val_score
from sklearn.metrics import r2_score
from scipy.stats import spearmanr, pearsonr
import matplotlib.pyplot as plt
from sklearn.model_selection import cross_val_score, cross_val_predict, ShuffleSplit
#from tpot import TPOTRegressor
import seaborn as sns
import random

# Load the data:
reactions = np.loadtxt("../data/fluxes.csv",delimiter=',', skiprows=1,usecols=0,dtype='str')
X_full = np.transpose(np.loadtxt("../data/transcripts.csv",delimiter=',', skiprows=1,usecols=range(1,53)))
y_full = np.transpose(np.loadtxt("../data/fluxes.csv",delimiter=',', skiprows=1,usecols=range(1,53)))
print(X_full.shape)

# Accuracy of a regressor is measured in terms of normalized root-mean-square error (NRMSE), with the normalization term used being the average flux. This accuracy term is also known as the coefficient of variation of the RMSD or CV(RMSD).  Lower values of NRMSE indicate less residual variance.

# Now loop over all measured metabolic fluxes which are greater than 200.0 umol/gdw/hr and create a random forest regressor for each flux:

samples = []
search = False
if search:
 for n in range(30):
  average_NRMSE = []
  samples_size = []
  for num_samples in range(4,52,4):
    flux_predicted_most_likely = []
    flux_predicted_most_lower = []
    flux_predicted_most_upper = []
    flux_predicted_NRMSE = []
    flux_measured =  []

    for i in range(2143):
#        y = np.copy(y_full[52-num_samples:52,i])
#        X = np.copy(X_full[52-num_samples:52,:])
        y = np.copy(y_full)
        X = np.copy(X_full)

        # Data is currently unshuffled; we should shuffle 
        # each X[i] with its corresponding y[i]
        perm = np.random.permutation(len(X))
        X = X[perm]
        y = y[perm]

        y = y[52-num_samples:52,i]
        X = X[52-num_samples:52,:]

#        print(num_samples,X.shape,y.shape)
#        print(num_samples,52-num_samples)
        if y.mean()>200.0:
            #use a random foreset regressor with 100 trees, using all 64 cores
            rf = RandomForestRegressor(n_estimators=100, oob_score=False, random_state=0,n_jobs=-1)

            #calculate cross validation using mean squared error in a 1:9 split
            scores = cross_val_score(rf,X,y=y,scoring='neg_mean_squared_error',cv=3)
            y_p = cross_val_predict(rf,X,y=y,cv=3)
            y_err = [ya - yp for ya,yp in zip(y,y_p)]
            
            #calculate normalized root-mean-square error
            NRMSE = np.sqrt(abs(scores.mean()))/y.mean()
            flux_predicted_NRMSE.append(NRMSE)
            
            #specify the training data in a 1:9 split of the 54 samples
            #test data set is the first 5 samples
            #training data set is the last 49 samples
            X_train = X[1:-1,:]
            y_train = y[1:-1]
            print(X_train.shape,'training dadta shape',)
            #Only look at strain number S. There are 54 total strains or scenarios.
            S = 0
            
            #make the random forest regressor
            rf.fit(X_train,y_train)
            print('Predicted flux for reaction',i,'in Strain',S,'is',rf.predict([X[S,:]])[0], 'while measured flux is',y[S])
            
            #grab the to be fuzzed transcripts for Strain S from X
            fuzzed_transcpits = np.array(X[S,:])

            possible_fluxes=[]
            
            #number of alternative measurements to search over that are
            #compatible with experimental error
            number_of_measurements = 1
            
            #assume experimental measurement error is 20%
            experimental_error = 0.20
            
            #first try max and min transcripts possible within experimental error
            fuzzed_transcpits_permuted = (1.0-experimental_error)*np.copy(fuzzed_transcpits)
            possible_fluxes.append(rf.predict(fuzzed_transcpits_permuted.reshape(-1, 7734))[0])
            fuzzed_transcpits_permuted = (1.0+experimental_error)*np.copy(fuzzed_transcpits)
            possible_fluxes.append(rf.predict(fuzzed_transcpits_permuted.reshape(-1, 7734))[0])
            
            #now loop over a large number of possible transcript measurements compatible with experimental error
            for j in range(number_of_measurements):
                fuzzed_transcpits_permuted = np.copy(fuzzed_transcpits)
                for k in range(len(fuzzed_transcpits_permuted)):
                    if random.choice([True, False]):
                        fuzzed_transcpits_permuted[k]= fuzzed_transcpits_permuted[k]*(1.0-experimental_error)
                    else:
                        fuzzed_transcpits_permuted[k]= fuzzed_transcpits_permuted[k]*(1.0+experimental_error)
                possible_fluxes.append(rf.predict(fuzzed_transcpits_permuted.reshape(-1, 7734))[0])
            possible_fluxes = np.array(possible_fluxes)
            
            flux_predicted_most_lower.append(possible_fluxes.min())
            flux_predicted_most_upper.append(possible_fluxes.max())
            #get predicted flux for Strain S
            flux_predicted_most_likely.append(rf.predict([X[S,:]])[0])
            
            flux_measured.append(y[S])
            '''            
            #plot predictions for flux i
            import seaborn as sns
            plt.scatter(y,y_p)
            plt.scatter(y,y)
            plt.title(str(num_samples)+' Samples: Reaction '+str(i))
            plt.xlabel('Measured flux [umol/gdw/hr]')
            plt.ylabel('Predicted flux [umol/gdw/hr]')
            plt.savefig(str(num_samples)+'Samples-Reaction'+str(i)+'.pdf')
            plt.close()
            plt.clf()
            '''
            print('done with reaction ',i,'. Number of samples is ',num_samples)
    '''
    my_range=range(1,len(flux_measured)+1)

    fig = plt.figure()
    ax = plt.subplot(111)

    # Shrink current axis by 20%
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])


    import seaborn as sns
    ax.plot(my_range, flux_predicted_most_likely,  "X",label='Most likely Predicted Flux')
    ax.plot(my_range, flux_predicted_most_lower,  "o",label='Lowest Predicted Flux')
    ax.plot(my_range, flux_predicted_most_upper,  "o",label='Highest Predicted Flux')
    ax.plot(my_range, flux_measured,  "X",label='Measured Flux')

    ax.set_title("Random Forest Transcripts to Fluxes Model", loc='left')
    ax.set_xlabel('Metabolic Reaction')
    ax.set_ylabel('Flux [umol/gdw/hr]')

    # Put a legend to the right of the current axis
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))

    #plt.show()
    plt.savefig(str(num_samples)+"Samples-predicted-fluxes.pdf")
    plt.cla()
    plt.clf()
    plt.close()

    my_range=range(1,len(flux_measured)+1)

    fig = plt.figure()
    ax = plt.subplot(111)

    # Shrink current axis by 20%
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])


    import seaborn as sns
    ax.plot(my_range, flux_predicted_NRMSE,  "X",label='Normalized Root-Mean-Square Error')

    ax.set_title("Random Forest Transcripts to Fluxes Model", loc='left')
    ax.set_xlabel('Metabolic Reaction')
    ax.set_ylabel('NRMSE (lower is better)')

    # Put a legend to the right of the current axis
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))

    #plt.show()
    plt.savefig(str(num_samples)+"Samples-nrmse.pdf")
    '''
    print('average nrmse is', np.average(flux_predicted_NRMSE),'for number of samples',num_samples)
    average_NRMSE.append(np.average(flux_predicted_NRMSE))
    samples_size.append(num_samples)


    '''  
    #plot histogram of error for biomass, i.e. the last trained regressor is on biomass
    sns.distplot(y_err)
    plt.title('Histogram of biomass residuals '+str(i))
    plt.xlabel('Flux error [umol/gdw/hr]')
    plt.ylabel('Frequency')
    plt.savefig(str(num_samples)+'Samples-histogram-biomass-residuals.pdf')
    plt.cla()
    plt.clf()
    plt.close()
    '''
  samples.append(average_NRMSE)
 np.savetxt('sample-size.csv',samples_size,delimiter=',',fmt='%i')
 np.savetxt('samples-NRMSE.csv',np.array(samples),delimiter=',')

df = pd.read_csv('samples-NRMSE.csv', header=None)
labels=(np.loadtxt('sample-size.csv',delimiter=',',dtype=str))

plt.cla()
plt.clf()
plt.close()
fig, ax = plt.subplots()
import seaborn as sns
#sns.pointplot(data=df)
sns.violinplot(data=df, inner=None,color='lightblue',alpha=.2)
sns.swarmplot(data=df, alpha=.5,color='black')



ax.set_xlabel('Number of samples for model training')
ax.set_ylabel('Average Accuracy [NRMSE]')
#plt.savefig("accuracy-vs-sample-size.pdf")
ax.set_xticklabels(labels)
plt.show()
