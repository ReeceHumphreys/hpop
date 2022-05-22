use finitediff::FiniteDiff;
use rgsl::legendre::associated_polynomials::legendre_sphPlm;
use std::fs::File;
use std::io::{prelude::*, BufReader, ErrorKind};
use std::string::String;

fn main() {
    let x = acceleration_gravity(6371000., 1.0482, 0.523699);
    let mag = f64::sqrt(x[0].powi(2) + x[1].powi(2) + x[2].powi(2));
    println!("{:?}", x);
    println!("{:?}", mag);
}

fn acceleration_gravity(r: f64, theta: f64, lambda: f64) -> Vec<f64> {
    const N_MAX: u32 = 50;
    const RADIUS_EARTH: f64 = 6378136.3; // [m]
    const GM_EARTH: f64 = 3986004.415e8; // [m^3/s^2]

    // r component
    let mut stuff = 0.;
    for n in 2..=N_MAX {
        let a = (n + 1) as f64 * (RADIUS_EARTH / r).powi(n as i32);
        let mut b = 0.;
        for m in 0..=n {
            let c = c_tesseral_coef(n, m);
            let s = s_tesseral_coef(n, m);
            let p_nm = legendre_sphPlm(n as i32, m as i32, theta.cos());
            b += ((c * f64::cos(m as f64 * lambda)) + (s * f64::sin(m as f64 * lambda))) * p_nm;
        }
        stuff = a * b;
    }
    let acceleration_r = (-GM_EARTH / r.powi(2)) * (1. + stuff);

    // theta component
    let mut stuff = 0.;
    for n in 2..=N_MAX {
        let a = (RADIUS_EARTH / r).powi(n as i32);
        let mut b = 0.;
        for m in 0..=n {
            let c = c_tesseral_coef(n, m);
            let s = s_tesseral_coef(n, m);
            let d_p_nm = |x: &Vec<f64>| -> f64 { legendre_sphPlm(n as i32, m as i32, x[1].cos()) };

            // Point at which gradient should be calculated
            let x: Vec<f64> = vec![r, theta, lambda];
            let grad_forward = x.central_diff(&d_p_nm);

            b += ((c * f64::cos(m as f64 * lambda)) + (s * f64::sin(m as f64 * lambda)))
                * grad_forward[1];
        }
        stuff = a * b;
    }
    let acceleration_theta = (GM_EARTH / r.powi(2)) * (stuff);

    // lambda component
    let mut stuff = 0.;
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
        stuff = a * b;
    }
    let acceleration_lambda = (GM_EARTH / (r.powi(2) * f64::sin(theta))) * (stuff);

    [acceleration_r, acceleration_theta, acceleration_lambda].to_vec()
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
