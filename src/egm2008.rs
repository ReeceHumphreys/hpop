use finitediff::FiniteDiff;
use ndarray::prelude::*;
use rgsl::legendre::associated_polynomials::legendre_sphPlm;
use std::fs::File;
use std::io::{prelude::*, BufReader, ErrorKind};
use std::string::String;

pub struct Model {
    file: File,
    n_max: usize,
}

impl Model {
    pub fn new(file_path: String, n_max: usize) -> Model {
        Model {
            file: File::open(file_path).expect("Could not find file"),
            n_max,
        }
    }

    pub fn acceleration_gravity(&mut self, r: f64, theta: f64, lambda: f64) -> Array1<f64> {
        let n_max: u32 = self.n_max as u32;
        const RADIUS_EARTH: f64 = 6378136.3; // [m]
        const GM_EARTH: f64 = 3986004.415e8; // [m^3/s^2]

        // r component
        let mut summation = 0.;
        for n in 2..=n_max {
            let a = (n + 1) as f64 * (RADIUS_EARTH / r).powi(n as i32);
            let mut b = 0.;
            for m in 0..=n {
                let c = self.c_tesseral_coef(n, m);
                let s = self.s_tesseral_coef(n, m);
                let p_nm = legendre_sphPlm(n as i32, m as i32, theta.cos());
                b += ((c * f64::cos(m as f64 * lambda)) + (s * f64::sin(m as f64 * lambda))) * p_nm;
            }
            summation += a * b;
        }
        let acceleration_r = (-GM_EARTH / r.powi(2)) * (1. + summation);

        // theta component
        let mut summation = 0.;
        for n in 2..=n_max {
            let a = (RADIUS_EARTH / r).powi(n as i32);
            let mut b = 0.;
            for m in 0..=n {
                let c = self.c_tesseral_coef(n, m);
                let s = self.s_tesseral_coef(n, m);

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
        for n in 2..=n_max {
            let a = (RADIUS_EARTH / r).powi(n as i32);
            let mut b = 0.;
            for m in 0..=n {
                let c = self.c_tesseral_coef(n, m);
                let s = self.s_tesseral_coef(n, m);
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

    /*
    -------------------------
    Parsing the EGM2008 data
    -------------------------
    */

    fn c_tesseral_coef(&mut self, l: u32, n: u32) -> f64 {
        let line: String = self.read_egm2008_file_line(l, n).unwrap();
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

    fn s_tesseral_coef(&mut self, l: u32, n: u32) -> f64 {
        let line: String = self.read_egm2008_file_line(l, n).unwrap();
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
    fn line_to_read(&mut self, l: u32, n: u32) -> Result<u32, std::io::Error> {
        if l < 2 {
            return Err(std::io::Error::new(
                ErrorKind::InvalidInput,
                "n must be greater than 1",
            ));
        }
        let x = 2..l;
        Ok(x.fold(l - 1, |a, b| a + b) + n - 1)
    }

    fn read_egm2008_file_line(&mut self, l: u32, n: u32) -> Result<String, std::io::Error> {
        let n: usize = self.line_to_read(l, n).unwrap() as usize;
        let buffer_reader = BufReader::new(&self.file);
        let result = buffer_reader.lines().nth(n);
        match result {
            Some(Ok(line)) => Ok(line),
            Some(Err(e)) => Err(e),
            None => Err(std::io::Error::new(ErrorKind::InvalidInput, "EOF reached")),
        }
    }
}
