mod egm2008;
use ndarray::prelude::*;
use ndarray_linalg::*;

fn main() {
    // Required data
    let rel_vel: Array1<f64> = array![1.0, 2.0, 3.0];
    let pos_sat: Array1<f64> = array![1.0, 2.0, 3.0];
    let pos_third: Array1<f64> = array![1.0, 2.0, 3.0];
    let pos_sat_to_sun: Array1<f64> = array![1.0, 2.0, 3.0];
    let mut egm_model = egm2008::Model::new("src/data/EGM2008_to2190_TideFree.txt".to_string(), 50);

    // Computing the acceleration from different effects.
    let a_grav = egm_model.acceleration_gravity(6371000. + (100. * 1e3), 1.0482, 0.523699);
    let a_drag = acceleration_drag(0.5, 0.5, 0.5, 0.5, &rel_vel);
    let a_third = acceleration_third_body(0.5, &pos_sat, &pos_third);
    let a_solar_rad =
        acceleration_solar_radiation_pressure(false, 0.1, 0.2, 0.3, 0.4, &pos_sat_to_sun);

    let _acceleration = &a_grav + a_drag + a_third + a_solar_rad;
    println!("{:?}", a_grav.norm_l2());
}

/// Returns the vector acceleration due to atmospheric drag.
///
/// # Arguments
///
/// * `density` - A f64 that is the density of the atmosphere.
/// * `drag_coefficient` - A f64 that is the drag coefficient of the object.
/// * `area` - A f64 that is the area of the object normal to the relative velocity.
/// * `mass` - A f64 that is the mass of the object.
/// * `relative_velocity` - A Vector3<f64> that is the relative velocity of the object.
fn acceleration_drag(
    density: f64,
    drag_coeff: f64,
    area: f64,
    mass: f64,
    relative_velocity: &Array1<f64>,
) -> Array1<f64> {
    let magnitude = -0.5 * density * (drag_coeff * area / mass) * relative_velocity.norm().powi(2);
    let direction = relative_velocity / relative_velocity.norm();
    direction * (magnitude)
}

/// Returns the vector acceleration due to a third body.
///
/// # Arguments
///
/// * `gravitational_parameter` - A f64 that is the gravitational parameter of 3rd body.
/// * `position_sat` - A Vector3<f64> that is the position vector of the object.
/// * `position_third_body` - A Vector3<f64> that is the position vector of the third body.
fn acceleration_third_body(
    gravitational_parameter: f64,
    position_sat: &Array1<f64>,
    position_third_body: &Array1<f64>,
) -> Array1<f64> {
    gravitational_parameter
        * ((position_sat / position_sat.norm().powi(3))
            - (position_third_body / position_third_body.norm().powi(3)))
}

/// Returns the vector acceleration due to a solar radiation pressure.
///
/// # Arguments
///
/// * `in_shadow` - A bool that indicates if the satellite is in shadow.
/// * `solar_radiation_coeff` - A f64 that is the coefficient of the solar radiation.
/// * `area` - A f64 that is the area of the object normal to the sun.
/// * `mass` - A f64 that is the mass of the object.
/// * `solar_pressure` - A f64 that is the solar pressure.
/// * `position_sat_to_sun` - A Vector3<f64> that is the position vector of the satellite to the sun.
fn acceleration_solar_radiation_pressure(
    in_shadow: bool,
    solar_radiation_coeff: f64,
    area: f64,
    mass: f64,
    solar_pressure: f64,
    position_sat_to_sun: &Array1<f64>,
) -> Array1<f64> {
    let in_shadow = in_shadow as i64;
    let scalar = (in_shadow as f64)
        * (solar_radiation_coeff * area / mass)
        * (solar_pressure / position_sat_to_sun.norm().powi(3));
    scalar * position_sat_to_sun
}
