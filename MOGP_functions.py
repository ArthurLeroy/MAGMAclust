import numpy as np
import mogptk

# create time array
n_points = 100
t = np.linspace(0.0, 6.0, n_points)

# channel 1
y1 = np.sin(6.0*t) + 0.2*np.random.normal(size=len(t))

# channel 2, phased version
y2 = np.sin(6.0*t + 2.0) + 0.2*np.random.normal(size=len(t))

# channel 3, added sinosoidal
y3 = np.sin(6.0*t) - np.sin(4.0*t) + 0.2*np.random.normal(size=len(t))

# channel 4, delayed and amplified
y4 = 3.0*np.sin(6.0 * (t-2.0)) + 0.3*np.random.normal(size=len(t))

# create dataset
dataset = mogptk.DataSet(
    mogptk.Data(t, y1, name='First channel'),
    mogptk.Data(t, y2, name='Second channel'),
    mogptk.Data(t, y3, name='Third channel'),
    mogptk.Data(t, y4, name='Fourth channel')
)

# remove 40% randomly
for data in dataset:
    data.remove_randomly(pct=0.4)

# remove second half of the first channel
dataset[0].remove_range(start=2.0)

# create model, uncomment for different kernels 
model = mogptk.MOSM(dataset, Q=2)
# model = mogptk.CSM(dataset, Q=2)
# model = mogptk.SM_LMC(dataset, Q=2)
# model = mogptk.CONV(dataset, Q=2)

# initialize parameters of kernel using LombScargle
model.init_parameters(method='LS', iters=500)

## Maybe use LBFGS for coherence
model.train(method='Adam', lr=0.1, iters=500,
  plot=True, error='MAE', verbose=True);

### From CSV ###

oil = mogptk.LoadCSV(filename='data/gonu/brent-daily.csv',
                     x_col='Date', y_col='Price', name='Oil')
 
gold = mogptk.LoadCSV(filename='data/gonu/lmba-gold-usd-am-daily.csv',
                      x_col='Date', y_col='Price', name='Gold',
                      na_values='.')

dataset = mogptk.DataSet(oil, gold)


for channel in dataset:
    # filter by date
    channel.filter('2015-01-01', '2018-12-31')
    
    # agregate from daily to weekly
    channel.aggregate('1W')
    
    # detrend the data set using a first-degree polynomial
    channel.transform(mogptk.TransformDetrend(degree=1))
    
model = mogptk.MOSM(dataset, Q=3)
model.init_parameters('SM')
model.train(method='Adam', lr=0.5, iters=500,
             verbose=True, plot=True, error='MAE');
             
model.plot_prediction();
    # remove by date range
    channel.remove_range('2016-11-15', '2017-01-01')
    channel.remove_randomly(pct=0.5)
    
    # set prediction range by date
    channel.set_prediction_range('2015-01-01', '2018-12-31', step='1D')
    
dataset.plot();

# save model
model.save('mosm_example')
# load model
loaded_model = mogptk.LoadModel('mosm_example')

loaded_model.plot_prediction(title='Loaded model');
