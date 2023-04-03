import numpy as np
import tensorflow as tf
import matplotlib.pyplot as plt

#Comment categories

# --> coding comments
'''
Machine learning theory comments

'''




# X-Data
N = 200
X = np.random.random(N)


sign = (- np.ones((N,)))**np.random.randint(2,size=N)
Y = np.sqrt(X) * sign




# Neural network: hiddens layers and ReLU activation function


act = tf.keras.layers.ReLU()

'''What are activation functions, and why are they necessary? 
Without an activation function like relu (also called a 
non-linearity), the Dense layer would consist of two linear 
operations—a dot product and an addition:
output = dot(W, input) + b
So the layer could only learn linear transformations (affine 
transformations) of the input data: the hypothesis space of 
the layer would be the set of all possible linear 
transformations of the input data into a 16-dimensional space. 
Such a hypothesis space is too restricted and wouldn’t benefit 
from multiple layers of representations, because a deep stack of 
linear layers would still implement a linear operation: adding 
more layers wouldn’t extend the hypothesis space.
In order to get access to a much richer hypothesis space that 
would benefit from deep representations, you need a non-linearity, 
or activation function. relu is the most popular activation 
function in deep learning, but there are many other 
candi- dates, which all come with similarly strange 
names: prelu, elu, and so on.

ReLU: Rectified Linear Unit activation
Very popular (most used circa 2017). Introduces non-linearity 
to model. Many different kinds of ReLU, their usefulness 
depending on the context. Avoided for a long time because 
it is not differentiable at x = 0. Researchers preferred 
"sigmoid" and "tanh". Formally it is a piecewise function 
connecting y = 0 and y = x at x = 0. It is easy to implement. 
It's only problem is that any negative value it takes in 
becomes zero immediately. Other ReLU have a way around this.

Why is it called "activation"? It's because at x = 0 the 
function is either activated or deactivated. It is just 
a function you use to get the node output, 
like "yes" or "no", or "ON" or "OFF".'''












nn_dp = tf.keras.models.Sequential([  # hidden layers here
	tf.keras.layers.Dense(10, activation=act, input_shape=(1,)), #activation function
	tf.keras.layers.Dense(10, activation=act),    
	tf.keras.layers.Dense(1, activation='linear')])

#A Sequential model is appropriate for a plain stack of layers where 
#each layer has exactly one input tensor and one output tensor.


'''
In neural networks, a hidden layer is located between the 
input and output of the algorithm, in which the function 
applies weights to the inputs and directs them through 
an activation function as the output. In short, the hidden 
layers perform nonlinear transformations of the inputs 
entered into the network. Hidden layers vary depending on 
the function of the neural network, and similarly, the 
layers may vary depending on their associated weights.



Hidden layers, simply put, are layers of mathematical 
functions each designed to produce an output specific 
to an intended result. For example, some forms of 
hidden layers are known as squashing functions.


'''






# Loss function
mse = tf.keras.losses.MeanSquaredError()
def loss_dp(y_true,y_pred):
 return mse(y_true,y_pred**2)

optimizer_dp = tf.keras.optimizers.Adam(learning_rate=0.001)
nn_dp.compile(optimizer=optimizer_dp, loss=loss_dp)


#Training
results_dp = nn_dp.fit(X, X, epochs=5, batch_size= 5, verbose=1)



# Plot 
plt.plot(X,Y,'.',label='Data points',color='lightgray')
plt.plot(X,nn_dp.predict(X),'.',label='Diff Phys',color='green')
plt.xlabel('y')
plt.ylabel('x')
plt.title('Differentiable Physics approach')
plt.legend()
plt.show()
