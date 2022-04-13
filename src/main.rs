use cgmath::prelude::*;
use cgmath::Vector3;
use finitediff::FiniteDiff;
use rgsl::legendre::associated_polynomials::legendre_sphPlm;
use std::fs::File;
use std::io::{prelude::*, BufReader, ErrorKind};
use std::string::String;

fn main() {
    let x = acceleration_gravity(63781369.3 + (1. * 1e3), 0.523699, 1.0482);
    //let mag = f64::sqrt(x[0].powi(2) + x[1].powi(2) + x[2].powi(2));
    println!("{:?}", x);
}

// /// Returns the vector acceleration due to atmospheric drag.
// ///
// /// # Arguments
// ///
// /// * `density` - A f64 that is the density of the atmosphere.
// /// * `drag_coefficient` - A f64 that is the drag coefficient of the object.
// /// * `area` - A f64 that is the area of the object normal to the relative velocity.
// /// * `mass` - A f64 that is the mass of the object.
// /// * `relative_velocity` - A Vector3<f64> that is the relative velocity of the object.
// fn acceleration_drag(
//     density: f64,
//     drag_coeff: f64,
//     area: f64,
//     mass: f64,
//     relative_velocity: Vector3<f64>,
// ) -> Vector3<f64> {
//     let magnitude = -0.5 * density * (drag_coeff * area / mass) * relative_velocity.magnitude2();
//     let direction = relative_velocity.normalize();
//     direction * (magnitude)
// }

// /// Returns the vector acceleration due to a third body.
// ///
// /// # Arguments
// ///
// /// * `gravitational_parameter` - A f64 that is the gravitational parameter of 3rd body.
// /// * `position_sat` - A Vector3<f64> that is the position vector of the object.
// /// * `position_third_body` - A Vector3<f64> that is the position vector of the third body.
// fn acceleration_third_body(
//     gravitational_parameter: f64,
//     position_sat: Vector3<f64>,
//     position_third_body: Vector3<f64>,
// ) -> Vector3<f64> {
//     gravitational_parameter
//         * ((position_sat / position_sat.magnitude().powi(3))
//             - (position_third_body / position_third_body.magnitude().powi(3)))
// }
fn delta_0(m: u64) -> u64 {
    if m == 0 {
        return 1;
    } else {
        return 0;
    }
}
fn normalization(n: u32, m: u32) -> f64 {
    let n = n as u64;
    let m = m as u64;
    let radicand: u64 = (factorial(n - m) * (2 * n + 1) * (2 - delta_0(m))) / factorial(n + m);
    f64::sqrt(radicand as f64)
}

fn factorial(num: u64) -> u64 {
    (1..=num).product()
}

fn acceleration_gravity(r: f64, theta: f64, lambda: f64) -> Vec<f64> {
    const N_MAX: u32 = 12;
    const RADIUS_EARTH: f64 = 6378136.3; // [m]
    const GM_EARTH: f64 = 3986004.415e8; // [m^3/s^2]

    // Define cost function `f(x)`
    // Gravitational potential of Earth using EGM2008
    let u = |x: &Vec<f64>| -> f64 {
        let r = x[0];
        let theta = x[1];
        let lambda = x[2];
        let mut b = 0_f64;
        for n in 2..N_MAX {
            let mut a: f64 = 0.;
            for m in 0..(n + 1) {
                let c = c_tesseral_coef(n, m);
                let s = s_tesseral_coef(n, m);
                let normalized_legendre =
                    normalization(n, m) * legendre_sphPlm(n as i32, m as i32, theta.cos());
                a += (c * f64::cos(m as f64 * lambda) + s * f64::sin(m as f64 * lambda))
                    * normalized_legendre;
            }
            b += (RADIUS_EARTH / r).powi(n as i32) * a;
        }
        let result = (GM_EARTH / r) * (1_f64 + b);
        println!("{:.64?}", result);
        result
    };

    // Point at which gradient should be calculated
    let x: Vec<f64> = vec![r, theta, lambda];

    // Calculate gradient of `u` at `x` using central differences
    let grad_forward = x.central_diff(&u);
    for point in &grad_forward {
        println!("{:.64?}", point);
    }
    grad_forward
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
// fn acceleration_solar_radiation_pressure(
//     in_shadow: bool,
//     solar_radiation_coeff: f64,
//     area: f64,
//     mass: f64,
//     solar_pressure: f64,
//     position_sat_to_sun: Vector3<f64>,
// ) -> Vector3<f64> {
//     let in_shadow = in_shadow as i64;
//     let scalar = (in_shadow as f64)
//         * (solar_radiation_coeff * area / mass)
//         * (solar_pressure / position_sat_to_sun.magnitude().powi(3));
//     scalar * position_sat_to_sun
// }

/*
-------------------------
Parsing the EGM2008 data
-------------------------
*/

// Is normalized already from data source
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
