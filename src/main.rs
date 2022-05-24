use finitediff::FiniteDiff;
use ndarray::prelude::*;
use ndarray_linalg::*;
use rgsl::legendre::associated_polynomials::legendre_sphPlm;
use std::fs::File;
use std::io::{prelude::*, BufReader, ErrorKind};
use std::string::String;

fn main() {
    // Required data
    let rel_vel: Array1<f64> = array![1.0, 2.0, 3.0];
    let pos_sat: Array1<f64> = array![1.0, 2.0, 3.0];
    let pos_third: Array1<f64> = array![1.0, 2.0, 3.0];
    let pos_sat_to_sun: Array1<f64> = array![1.0, 2.0, 3.0];

    // Computing the acceleration from different effects.
    let a_grav = acceleration_gravity(6371000. + (100. * 1e3), 1.0482, 0.523699);
    let a_drag = acceleration_drag(0.5, 0.5, 0.5, 0.5, &rel_vel);
    let a_third = acceleration_third_body(0.5, &pos_sat, &pos_third);
    let a_solar_rad =
        acceleration_solar_radiation_pressure(false, 0.1, 0.2, 0.3, 0.4, &pos_sat_to_sun);

    let acceleration = a_grav + a_drag + a_third + a_solar_rad;
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

fn acceleration_gravity(r: f64, theta: f64, lambda: f64) -> Array1<f64> {
    const N_MAX: u32 = 50;
    const RADIUS_EARTH: f64 = 6378136.3; // [m]
    const GM_EARTH: f64 = 3986004.415e8; // [m^3/s^2]

    // r component
    let mut summation = 0.;
    for n in 2..=N_MAX {
        let a = (n + 1) as f64 * (RADIUS_EARTH / r).powi(n as i32);
        let mut b = 0.;
        for m in 0..=n {
            let c = c_tesseral_coef(n, m);
            let s = s_tesseral_coef(n, m);
            let p_nm = legendre_sphPlm(n as i32, m as i32, theta.cos());
            b += ((c * f64::cos(m as f64 * lambda)) + (s * f64::sin(m as f64 * lambda))) * p_nm;
        }
        summation += a * b;
    }
    let acceleration_r = (-GM_EARTH / r.powi(2)) * (1. + summation);

    // theta component
    let mut summation = 0.;
    for n in 2..=N_MAX {
        let a = (RADIUS_EARTH / r).powi(n as i32);
        let mut b = 0.;
        for m in 0..=n {
            let c = c_tesseral_coef(n, m);
            let s = s_tesseral_coef(n, m);

            // The partial derivative of the legendre polynomial with respect to theta
            let u = |x: &Vec<f64>| -> f64 { legendre_sphPlm(n as i32, m as i32, x[1].cos()) };
            let x: Vec<f64> = vec![r, theta, lambda];
            let d_pnm_d_theta = x.central_diff(&u)[1];

            b += ((c * f64::cos(m as f64 * lambda)) + (s * f64::sin(m as f64 * lambda)))
                * d_pnm_d_theta;
        }
        summation += a * b;
    }
    let acceleration_theta = (GM_EARTH / r.powi(2)) * (summation);

    // lambda component
    let mut summation = 0.;
    for n in 2..=N_MAX {
        let a = (RADIUS_EARTH / r).powi(n as i32);
        let mut b = 0.;
        for m in 0..=n {
            let c = c_tesseral_coef(n, m);
            let s = s_tesseral_coef(n, m);
            let p_nm = legendre_sphPlm(n as i32, m as i32, theta.cos());
            b += m as f64
                * ((-c * f64::sin(m as f64 * lambda)) + (s * f64::cos(m as f64 * lambda)))
                * p_nm;
        }
        summation += a * b;
    }
    let acceleration_lambda = (GM_EARTH / (r.powi(2) * f64::sin(theta))) * (summation);

    array![acceleration_r, acceleration_theta, acceleration_lambda]
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

/*
-------------------------
Parsing the EGM2008 data
-------------------------
*/

fn c_tesseral_coef(l: u32, n: u32) -> f64 {
    let line: String = read_egm2008_file_line(l, n).unwrap();
    let terms: Vec<&str> = line.split_whitespace().collect();
    let c_term: f64 = terms[2]
        .to_string()
        .chars()
        .map(|x| match x {
            'D' => 'e',
            _ => x,
        })
        .collect::<String>()
        .parse::<f64>()
        .unwrap();
    c_term
}

fn s_tesseral_coef(l: u32, n: u32) -> f64 {
    let line: String = read_egm2008_file_line(l, n).unwrap();
    let terms: Vec<&str> = line.split_whitespace().collect();
    let s_term: f64 = terms[3]
        .to_string()
        .chars()
        .map(|x| match x {
            'D' => 'e',
            _ => x,
        })
        .collect::<String>()
        .parse::<f64>()
        .unwrap();
    s_term
}

// Outputs the line of the EGM2008 file to rad for a given n and m.
fn line_to_read(l: u32, n: u32) -> Result<u32, std::io::Error> {
    if l < 2 {
        return Err(std::io::Error::new(
            ErrorKind::InvalidInput,
            "n must be greater than 1",
        ));
    }
    let x = 2..l;
    Ok(x.fold(l - 1, |a, b| a + b) + n - 1)
}

fn read_egm2008_file_line(l: u32, n: u32) -> Result<String, std::io::Error> {
    let file = File::open("src/data/EGM2008_to2190_TideFree.txt").unwrap();
    let buf_reader = BufReader::new(file);
    let line_to_read: usize = line_to_read(l, n).unwrap() as usize;
    let u = buf_reader
        .lines()
        .nth(line_to_read)
        .expect("Could not read line")
        .unwrap();
    Ok(u)
}
