import numpy as np 
import pandas as pd

def load_data(file_path, circle_points, n_slices):

    #load data from csv file
    df = pd.read_csv(file_path) 

    #checks if file is in the correct format
    total_expected_rows = circle_points * n_slices 

    if len(df) != total_expected_rows:
        raise ValueError(f"Expected {total_expected_rows} rows, but got {len(df)}. Please check the input file.")
    
    #make a copy of the dataframe to work with
    df_worked = df.copy()

    #centrepoints and wall coordinates for each slice
    wall_points = df_worked[['X', 'Y', 'Z']].values
    centrepoints = df_worked[['X_centreline','Y_centreline','Z_centreline']].values

    #the centre points are repeated for each slice, so we take only the first occurrence for each slice
    unique_centrepoints = df_worked.iloc[::circle_points][['X_centreline','Y_centreline','Z_centreline']].values

    return df_worked, wall_points, unique_centrepoints, centrepoints


def axial_vectors(unique_centrepoints):
    
    #calculating the differences between consecutive centre points to get the axial vectors direction
    axial_vectors = np.diff(unique_centrepoints, axis=0)

    norms = np.linalg.norm(axial_vectors, axis=1, keepdims=True)
    norms[norms == 0] = 1  # normalising the vectors to get unit vectors
    axial_unit_vectors = axial_vectors / norms

    #we need to calculate the axial vector for the last slice as well, which is the same as the second last one
    last_vector = axial_unit_vectors[-1]
    axial_unit_vectors = np.vstack([axial_unit_vectors, last_vector])

    return axial_unit_vectors


def radial_vectors(wall_points,centrepoints):

    radial_vec = wall_points - centrepoints
    norms = np.linalg.norm(radial_vec, axis=1, keepdims=True)
    norms[norms == 0] = 1  # normalising the vectors to get unit vectors

    return radial_vec / norms


def tangential_vectors(axial_vectors, radial_vectors):
    # Compute tangential vectors as the cross product of axial and radial vectors
    return np.cross(axial_vectors, radial_vectors)


def process_data(file_path, circle_points, n_slices, output_file = None):
   
    print(f"reading file: {file_path}...")

    #loading data from csv file
    df, wall_points, centrepoints, unique_centrepoints = load_data(file_path, circle_points, n_slices)
    print("calculating vectors...")

    #calculating axial, radial and tangential vectors
    axial_vector_slice = axial_vectors(unique_centrepoints)
    axial_vectors = np.repeat(axial_vector_slice, circle_points, axis=0)

    radial_vec = radial_vectors(wall_points, centrepoints)
    tangential_vec = tangential_vectors(axial_vectors, radial_vec)

    df['nv_x'] = radial_vec[:, 0]
    df['nv_y'] = radial_vec[:, 1]
    df['nv_z'] = radial_vec[:, 2]   

    df['tv_x'] = tangential_vec[:, 0]
    df['tv_y'] = tangential_vec[:, 1]
    df['tv_z'] = tangential_vec[:, 2]
    
    if output_file:
        df.to_csv(output_file, index=False)
        print(f"processed data saved to: {output_file}")

    return df



