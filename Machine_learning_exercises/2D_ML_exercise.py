
from mpl_toolkits import mplot3d
import numpy as np
import matplotlib.pyplot as plt
import numpy as np
import tensorflow as tf
from matplotlib import cm



# to get parabola to fit better try two things
# 1: Nonuniform sampling: Get more points by y=0
# 2: Change NN training (vague)



#Thoughts
# Try another nonlinear activation function

fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
N = 100



sign = (- np.ones((N,)))**np.random.randint(2,size=N)
X = np.zeros(N,float);Y = np.zeros(N,float)

for i in range(0,len(X)):
 X[i] = np.random.random()
 Y[i] = np.random.random()


R = np.sqrt(X) + np.sqrt(Y)
Z = sign*R

print(Z.shape,R.shape)


# A StrMethodFormatter is used automatically

# Add a color bar which maps values to colors.


act = tf.keras.layers.ReLU()
nn_dp = tf.keras.models.Sequential([
	tf.keras.layers.Dense(10, activation=act, input_shape=(1,)),
	tf.keras.layers.Dense(10, activation=act),
	tf.keras.layers.Dense(10, activation=act),
 	tf.keras.layers.Dense(10, activation=act),
	tf.keras.layers.Dense(10, activation=act),
	tf.keras.layers.Dense(1, activation='linear')])


# Loss function
mse = tf.keras.losses.MeanSquaredError()
def loss_dp(z_true,z_pred):
 return mse(z_true,z_pred**2)

optimizer_dp = tf.keras.optimizers.Adam(learning_rate=0.001)
nn_dp.compile(optimizer=optimizer_dp, loss=loss_dp)



#Training
results_dpx = nn_dp.fit(R,Z, epochs=100, batch_size= 5, verbose=1)
#results_dpy = nn_dp.fit(Y, Y, epochs=10, batch_size= 5, verbose=1)
z = nn_dp.predict(Z)


# Plot solution









surf1 = ax.scatter(X, Y, Z,alpha=0.3)

surf2 = ax.scatter(X,Y,z)

plt.title('Differentiable Physics approach')
plt.show()
