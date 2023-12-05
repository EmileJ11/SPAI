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
rir_data = scipy.io.loadmat('C:/Users/jnema/mat_files/speakers_all_discrete_positions.mat')['speakers_all_discrete_positions']


# print("DOA data shape:", doa_data.shape)
# print("Pressure data shape:", pressure_data.shape)
# print("RIR data shape:", rir_data.shape)

# print("DOA data type:", doa_data.dtype)
# print("Pressure data type:", pressure_data.dtype)
# print("RIR data type:", rir_data.dtype)

# # Optionally, print the first few rows to visually inspect the data
# print("DOA data sample:", doa_data[:5])
# print("Pressure data sample:", pressure_data[:5])
# print("RIR data sample:", rir_data[:5])

def flatten_nested_arrays(data):
    # Assuming each element in the data is an array of the same shape
    flattened = np.concatenate(data.ravel())
    return flattened

# Flatten DOA data
doa_data_flattened = flatten_nested_arrays(doa_data)

# Flatten Pressure data
pressure_data_flattened = flatten_nested_arrays(pressure_data)

# Flattem RIR data
rir_data_flattened = flatten_nested_arrays(rir_data)

# # Ensure the flattened data is in the right shape
# print("Flattened DOA data shape:", doa_data_flattened.shape)
# print("Flattened Pressure data shape:", pressure_data_flattened.shape)




# Initialize your PINN model
model = PINN(num_neurons_per_layer=20, num_outputs=1) 


# Normalize DOA data
scaler_doa = StandardScaler()
doa_data_scaled = scaler_doa.fit_transform(doa_data_flattened)

# Normalize Pressure data
scaler_pressure = StandardScaler()
pressure_data_scaled = scaler_pressure.fit_transform(pressure_data_flattened)

# Concatenate the scaled DOA and Pressure data
input_data_scaled = np.concatenate([doa_data_scaled, pressure_data_scaled], axis=1)

# Normalize RIR data
scaler_rir = StandardScaler()
rir_data_scaled = scaler_rir.fit_transform(rir_data.reshape(-1, 1) if rir_data.ndim == 1 else rir_data_flattened)




#Split the data into training and validation sets
x_train, x_val, y_train, y_val = train_test_split(input_data_scaled, rir_data_scaled, test_size=0.2, random_state=42)


def pde_residuals(model, collocation_points_np, c):
    """
    Calculate the residuals of the wave equation (PDE) at collocation points.

    Parameters:
    model: The neural network model (PINN).
    collocation_points_np: Points where the PDE is evaluated (time and spatial coordinates) as a NumPy array.
    c: Speed of sound.

    Returns:
    Residuals of the PDE.
    """
    # Convert collocation points to a TensorFlow tensor
    collocation_points = tf.convert_to_tensor(collocation_points_np, dtype=tf.float32)

    with tf.GradientTape(persistent=True) as tape:
        tape.watch(collocation_points)

        # Predict the response using the neural network
        p = model(collocation_points)

        # Extract time and spatial coordinates from collocation points
        t = collocation_points[:, 0]


        # Check if p depends on t
        if tape.gradient(p, t) is None:
            print("p does not depend on t. Check model architecture and input.")
            return 0

        # Compute first and second time derivatives
        dp_dt = tape.gradient(p, t)
        if dp_dt is None:
            print("dp_dt is None. Gradient computation failed.")
            return 0

        d2p_dt2 = tape.gradient(dp_dt, t)
        if d2p_dt2 is None:
            print("d2p_dt2 is None. Gradient computation failed.")
            return 0


        x = collocation_points[:, 1]
        y = collocation_points[:, 2]
        z = collocation_points[:, 3]

        # Compute first and second time derivatives
        dp_dt = tape.gradient(p, t)
        d2p_dt2 = tape.gradient(dp_dt, t)

        # Compute spatial derivatives (Laplacian)
        dp_dx = tape.gradient(p, x)
        d2p_dx2 = tape.gradient(dp_dx, x)
        dp_dy = tape.gradient(p, y)
        d2p_dy2 = tape.gradient(dp_dy, y)
        dp_dz = tape.gradient(p, z)
        d2p_dz2 = tape.gradient(dp_dz, z)

        laplacian_p = d2p_dx2 + d2p_dy2 + d2p_dz2

    # Residual of the wave equation
    wave_eq_residual = d2p_dt2 - (1 / c**2) * laplacian_p

    return tf.reduce_mean(tf.square(wave_eq_residual))


def custom_loss_function(y_true, y_pred, model, collocation_points, c, epsilon_f, epsilon_d):
    """
    Custom loss function for the PINN model, combining Ldata and LPDE.

    Parameters:
    y_true: The actual target values from your data (RIR data).
    y_pred: The predicted values from the model.
    model: The PINN model.
    collocation_points: Points to evaluate the PDE.
    c: Speed of sound.
    epsilon_f: Adaptive weight for the PDE loss.
    epsilon_d: Adaptive weight for the data fidelity loss.

    Returns:
    Total loss value.
    """
    # Data fidelity loss (Ldata)
    data_fidelity_loss = tf.reduce_mean(tf.square(y_true - y_pred))

    # PDE loss (LPDE) based on the wave equation
    pde_loss_value = pde_residuals(model, collocation_points, c)

    # Combine losses with adaptive weights
    total_loss = (1 / (2 * epsilon_f**2)) * pde_loss_value + \
                 (1 / (2 * epsilon_d**2)) * data_fidelity_loss + \
                 tf.math.log(epsilon_f * epsilon_d)

    return total_loss

# Initialize adaptive weights
epsilon_f = 1.0  # Initial value for physics loss weight
epsilon_d = 1.0  # Initial value for data fidelity loss weight
speed_of_sound = 343


# Extract min and max values for x, y, z from doa_data
x_values = np.hstack(doa_data_flattened[:, 0])  # Flatten any nested structures in the first column
y_values = np.hstack(doa_data_flattened[:, 1])  # Flatten any nested structures in the second column
z_values = np.hstack(doa_data_flattened[:, 2])  # Flatten any nested structures in the third column

x_min, x_max = np.min(x_values), np.max(x_values)
y_min, y_max = np.min(y_values), np.max(y_values)
z_min, z_max = np.min(z_values), np.max(z_values)



t_min, t_max = np.min(rir_data_flattened), np.max(rir_data_flattened)

# print('X min: ', x_min)
# print('T min: ', t_min)
# print('T max: ', t_max)


# Define the range for space (x, y, z) and time (t)
# Replace these with the actual ranges for your problem
x_range = (x_min, x_max)
y_range = (y_min, y_max)
z_range = (z_min, z_max)
t_range = (t_min, t_max)

# Number of points in each dimension
num_points_space = int(np.cbrt(12000))  # Cube root of 12000 for even distribution in 3D space
num_points_time = num_points_space

# Generate uniform points in each dimension
x = np.linspace(x_range[0], x_range[1], num_points_space)
y = np.linspace(y_range[0], y_range[1], num_points_space)
z = np.linspace(z_range[0], z_range[1], num_points_space)
t = np.linspace(t_range[0], t_range[1], num_points_time)

# Create a grid of points (collocation points)
xx, yy, zz, tt = np.meshgrid(x, y, z, t)
collocation_points = np.array([xx.ravel(), yy.ravel(), zz.ravel(), tt.ravel()]).T



# Compile the model with the custom loss function
model.compile(optimizer='adam', loss=lambda y_true, y_pred: custom_loss_function(y_true, y_pred, model, collocation_points, speed_of_sound, epsilon_f, epsilon_d))

# Train the model
history = model.fit(x_train, y_train, epochs=1, batch_size=32, validation_data=(x_val, y_val))  # Adjust epochs and batch_size as needed


predicted_y = model.predict(x_val)


predicted_y = predicted_y.squeeze()  # Remove extra dimensions if necessary
plt.scatter(y_val, predicted_y)
plt.xlabel('Actual')
plt.ylabel('Predicted')
plt.title('Predictions vs Actual Data')

# plt.subplot(1, 2, 2)
# plt.plot(history.history['loss'], label='Training loss')
# plt.xlabel('Epoch')
# plt.ylabel('Loss')
# plt.title('Training Loss')
# plt.legend()
# plt.show()
