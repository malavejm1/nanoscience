import numpy as np
import tensorflow as tf
import matplotlib.pyplot as plt

# to get parabola to fit better try two things
# 1: Nonuniform sampling: Get more points by y=0
# 2: Change NN training (vague)



#Thoughts
# Try another nonlinear activation function

fig=plt.figure()


N = 200
X1 = np.random.random(N)
X2 = np.linspace(0,0.1,num=N)
X = np.zeros(len(X1)+len(X2),float)
for i in range(0,len(X1)):
 X[i] = X1[i]
 X[i+len(X2)] = X2[i]



sign = (- np.ones((len(X),)))**np.random.randint(2,size=len(X))
Y = np.sqrt(X) * sign
Z = np.sqrt(X) * sign




act = tf.keras.layers.ReLU()
nn_dp = tf.keras.models.Sequential([
	tf.keras.layers.Dense(10, activation=act, input_shape=(1,)),
	tf.keras.layers.Dense(10, activation=act),
	tf.keras.layers.Dense(10, activation=act),
	tf.keras.layers.Dense(1, activation='linear')])


# Loss function
mse = tf.keras.losses.MeanSquaredError()
def loss_dp(y_true,y_pred):
 return mse(y_true,y_pred**2)

optimizer_dp = tf.keras.optimizers.Adam(learning_rate=0.001)
nn_dp.compile(optimizer=optimizer_dp, loss=loss_dp)


#Training
results_dp = nn_dp.fit(X, X, epochs=10, batch_size= 5, verbose=1)



# Plot shitty solution
plt.plot(X,Y,'.',label='Data points',color='lightgray')
plt.plot(X,nn_dp.predict(X),'.',label='Diff Phys',color='green')
plt.xlabel('y')
plt.ylabel('x')
plt.title('Differentiable Physics approach')
plt.legend()
plt.show()
