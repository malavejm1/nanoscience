import numpy as np
import tensorflow as tf
import matplotlib.pyplot as plt

# X-Data
N = 200 
X = np.random.random(N)




# Geneation Y-data
sign = (- np.ones((N,)))**np.random.randint(2,size=N) 
Y = np.sqrt(X) * sign

# Neural network: three hidden layers and ReLU activations
act = tf.keras.layers.ReLU()
nn_sv = tf.keras.models.Sequential([
	tf.keras.layers.Dense(10, activation=act, input_shape=(1,)),
	tf.keras.layers.Dense(10, activation=act),
	tf.keras.layers.Dense(1, activation='linear')])


# Loss function
loss_sv = tf.keras.losses.MeanSquaredError()
optimizer_sv = tf.keras.optimizers.Adam(learning_rate=0.001)
nn_sv.compile(optimizer=optimizer_sv, loss=loss_sv)

#Training
results_sv = nn_sv.fit(X, Y, epochs=5, batch_size= 5, verbose=1)



# Plot shitty solution
plt.plot(X,Y,'.',label='Data points',color='lightgray')
plt.plot(X,nn_sv.predict(X),'.',label='Supervised',color='red')
plt.xlabel('y')
plt.ylabel('x')
plt.title('Standard approach')
plt.legend()
plt.show()
