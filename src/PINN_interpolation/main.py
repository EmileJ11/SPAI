import tensorflow as tf
import numpy as np
import matplotlib.pyplot as plt
import scipy.io
from pinn_model import PINN
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler


# Load the data
doa_data = scipy.io.loadmat('C:/Users/jnema/mat_files/DOA.mat')['DOA']
pressure_data = scipy.io.loadmat('C:/Users/jnema/mat_files/pressure.mat')['P']


# Extracting DOA data for speaker 1 at points 1, 2, and 3
doa_point_1_speaker_1 = doa_data[0][0]  # First array of the first pair
doa_point_2_speaker_1 = doa_data[1][0]  
doa_point_3_speaker_1 = doa_data[2][0]  

pressure_point_1_speaker_1 = pressure_data[0][0]
pressure_point_2_speaker_1 = pressure_data[1][0]
pressure_point_3_speaker_1 = pressure_data[2][0]



# To test:
# print("Shape of P1 S1: ", doa_point_1_speaker_1.shape)
# print("Shape of DOA: ", doa_data.shape)


# Concatenate the scaled DOA and Pressure data
combined_point_1 = np.concatenate([doa_point_1_speaker_1, pressure_point_1_speaker_1], axis=1)
combined_point_2 = np.concatenate([doa_point_2_speaker_1, pressure_point_2_speaker_1], axis=1)
combined_point_3 = np.concatenate([doa_point_3_speaker_1, pressure_point_3_speaker_1], axis=1)


sampling_rate = 48000  # 48000 samples per second
speed_of_sound = 343
time_interval = 1 / sampling_rate  # Time interval between samples
num_samples = pressure_data.shape[0]  # Number of rows in your data
row_index = np.arange(num_samples)
time_index = row_index / sampling_rate



def custom_loss_function(y_true, y_pred, model, collocation_tensor, epsilon_f, epsilon_d):
    # Data fidelity loss (e.g., Mean Squared Error)
    data_fidelity_loss = tf.reduce_mean(tf.square(y_true - y_pred))

    # Physics-informed term (e.g., wave equation residual)
    # Implement the function to calculate this term
    physics_loss = calculate_physics_loss(model, collocation_tensor, speed_of_sound)

    # Combine the losses with adaptive weighting
    total_loss = (1 / (2 * epsilon_f**2)) * physics_loss + \
                 (1 / (2 * epsilon_d**2)) * data_fidelity_loss + \
                 tf.math.log(epsilon_f * epsilon_d)

    return total_loss


def calculate_physics_loss(model, collocation_tensor, c):

    with tf.GradientTape(persistent=True) as tape:
        tape.watch(collocation_tensor)
        predictions = model(collocation_tensor)

        # Split the predictions into DOA and pressure components
        predicted_doa = predictions[:, :3]  # Assuming first 3 columns are DOA
        predicted_pressure = predictions[:, 3]  # Assuming last column is pressure

        # Calculate the wave equation residual for pressure
        pressure_residual = calculate_pressure_residual(predicted_pressure, collocation_tensor, c)

        # Add constraints for DOA consistency or other physical laws governing DOA
        doa_residual = calculate_doa_residual(predicted_doa)

    # Combine residuals for total physics loss
    return tf.reduce_mean(tf.square(pressure_residual)) + tf.reduce_mean(tf.square(doa_residual))


def find_min_max_coordinates(doa_points):
    # Concatenate the DOA points
    all_points = np.concatenate([doa_points[0], doa_points[1], doa_points[2]], axis=0)

    # Find min and max for x, y, and z
    x_min, x_max = np.min(all_points[:, 0]), np.max(all_points[:, 0])
    y_min, y_max = np.min(all_points[:, 1]), np.max(all_points[:, 1])
    z_min, z_max = np.min(all_points[:, 2]), np.max(all_points[:, 2])

    return x_min, x_max, y_min, y_max, z_min, z_max



x_min, x_max, y_min, y_max, z_min, z_max = find_min_max_coordinates([doa_point_1_speaker_1, doa_point_2_speaker_1, doa_point_3_speaker_1])



def generate_collocation_points():
    # Define the range of x, y, z, and t
    x_range = (x_min, x_max)
    y_range = (y_min, y_max)
    z_range = (z_min, z_max)
    t_range = (0, 450000)

    # Sample points within these ranges
    x_points = np.linspace(*x_range, num=50)
    y_points = np.linspace(*y_range, num=50)
    z_points = np.linspace(*z_range, num=50)
    t_points = np.linspace(*t_range, num=50)

    # Create a grid of points
    xx, yy, zz, tt = np.meshgrid(x_points, y_points, z_points, t_points)

    collocation_tensor = tf.convert_to_tensor(np.array([xx.ravel(), yy.ravel(), zz.ravel(), tt.ravel()]).T, dtype=tf.float32)

    return collocation_tensor



def calculate_pressure_residual(pressure, collocation_points, c):
    with tf.GradientTape(persistent=True) as tape:
        time_tensor = tf.convert_to_tensor(time_index, dtype=tf.float32)
        tape.watch(time_tensor)

        pressure_tensor = tf.convert_to_tensor(pressure, dtype=tf.float32)
        tape.watch(pressure_tensor)

        # Compute first and second time derivatives
        dp_dt = tape.gradient(pressure_tensor, time_tensor)
        d2p_dt2 = tape.gradient(dp_dt, time_tensor)

        dp_dx = tape.gradient(pressure, collocation_points[:, 0])
        d2p_dx2 = tape.gradient(dp_dx, collocation_points[:, 0])
        dp_dy = tape.gradient(pressure, collocation_points[:, 1])
        d2p_dy2 = tape.gradient(dp_dy, collocation_points[:, 1])
        dp_dz = tape.gradient(pressure, collocation_points[:, 2])
        d2p_dz2 = tape.gradient(dp_dz, collocation_points[:, 2])

        laplacian_p = d2p_dx2 + d2p_dy2 + d2p_dz2

    wave_eq_residual = d2p_dt2 - (speed_of_sound**2) * laplacian_p

    return wave_eq_residual


def calculate_doa_residual(doa):
    normalized_doa = tf.nn.l2_normalize(doa, axis=1)
    residual = tf.reduce_mean(tf.square(doa - normalized_doa))

    return residual


def train_model(X_train, y_train):
    # Initialize the PINN model
    model = PINN(num_neurons_per_layer=20, num_outputs=4)

    # Compile the model with an optimizer and the custom loss function
    model.compile(optimizer='adam', 
                  loss=lambda y_true, y_pred: custom_loss_function(y_true, y_pred, model, generate_collocation_points(), 1.0, 1.0))

    # Train the model
    model.fit(X_train, y_train, epochs=100, batch_size=32, validation_split=0.2)

    return model





# Extract the first 1000 rows for training and testing
train_data = np.concatenate([combined_point_1[:1000, :], combined_point_3[:1000, :]])
test_data = combined_point_2[:1000, :]

# Splitting the train data into input (X) and output (y)
X_train = train_data[:, :4]  # First 4 columns: x, y, z, and pressure
y_train = train_data[:, :4]  # Output is also x, y, z, and pressure

# Testing data
X_test = test_data[:, :4]
y_test = test_data[:, :4]

# Training the model
# Note: You need to define the 'train_model' function based on your training procedure
model = train_model(X_train, y_train)

# Predicting
predictions = model.predict(X_test)

# Calculating average error for DOA and pressure
average_error_doa = np.mean(np.abs(y_test[:, :3] - predictions[:, :3]), axis=0)  # For DOA x, y, z
average_error_pressure = np.mean(np.abs(y_test[:, 3] - predictions[:, 3]))  # For pressure

print("Average Error in DOA (x, y, z):", average_error_doa)
print("Average Error in Pressure:", average_error_pressure)
