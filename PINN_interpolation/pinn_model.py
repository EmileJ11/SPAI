import tensorflow as tf

class PINN(tf.keras.Model):
    def __init__(self, num_neurons_per_layer, num_outputs):
        super(PINN, self).__init__()
        self.dense1 = tf.keras.layers.Dense(num_neurons_per_layer, activation=tf.sin)
        self.dense2 = tf.keras.layers.Dense(num_neurons_per_layer, activation=tf.sin)
        self.out = tf.keras.layers.Dense(num_outputs, activation=None)  # Output layer

    def call(self, inputs):
        x = self.dense1(inputs)
        x = self.dense2(x)
        return self.out(x)