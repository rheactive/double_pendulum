use macroquad::prelude::*;

use std::{thread, time};
use std::fs;
use std::io::{Write};
use std::path::{PathBuf};

const PI: f64 = 3.1415926536;
const SPEED: f64 = 1.0;

// express angle in degrees
fn angle_in_degrees (angle: f64) -> f64 {
    let in_degrees;
    if angle > 2.0 * PI {
        let full_angles = (angle / (2.0 * PI)) as u32;
        in_degrees = angle * 180.0 / PI - ((360 * full_angles) as f64)
    } else {
        if angle < - 2.0 * PI {
            let full_angles = (- angle / (2.0 * PI)) as u32;
            in_degrees = angle * 180.0 / PI + ((360 * full_angles) as f64)
        } else {
            in_degrees = angle * 180.0 / PI
        }
    }
    if in_degrees < 0.0 {360.0 + in_degrees}
    else {in_degrees}
}

// function for the right hand side of the equation
fn right_hand_side (theta: f64, phi: f64, p: f64, q: f64, omega2: f64) -> [f64; 5] {
    let cs = (theta - phi).cos();
    let sn = (theta - phi).sin();
    let a = 3.0 / (16.0 - 9.0 * cs * cs);
    let theta_dot = a * (8.0 * p - 3.0 * q * cs);
    let phi_dot = a * (2.0 * q - 3.0 * p * cs);
    let b = theta_dot * phi_dot * sn;
    let p_dot = - b - 3.0 * omega2 * theta.sin();
    let q_dot = b - omega2 * phi.sin();
    // Also evaluate total energy
    let h = p * theta_dot + q * phi_dot - 2.0 * omega2 * (3.0 * theta.cos() + phi.cos());
    [theta_dot, phi_dot, p_dot, q_dot, h]
}

fn get_window_par (height: f32, width: f32) -> [f32; 4] {
    let l = height / 5.0;
    let w = width / 100.0;
    let x0 = width / 2.0;
    let y0 = height / 2.0;
    [l, w, x0, y0]
}

fn draw_pendulum (theta: f64, phi: f64, par: [f32; 4]) {
    // calculate and plot pendulum positions
    let l = par[0];
    let w = par[1];
    let x0 = par[2];
    let y0 = par[3];

    let x1 = x0 - l * (theta.sin() as f32);
    let y1 = y0 + l * (theta.cos() as f32);
    let x2 = x1 - l * (phi.sin() as f32);
    let y2 = y1 + l * (phi.cos() as f32);

    draw_line(x0, y0, x1, y1, w, BLACK);
    draw_line(x1, y1, x2, y2, w, BLACK);
    draw_circle(x0, y0, w, BLACK);
    draw_circle(x1, y1, w, BLUE);
    draw_circle(x2, y2, w, RED);
}

#[macroquad::main("Double Pendulum")]
async fn main() {

    //make directory
    let dir_name = "results";
    fs::create_dir_all(dir_name).expect("Error creating directory");
	
    // file for data
    //file for nanostructures data
    let fl_name = "double_pendulum.dat";
        let file_path: PathBuf = [dir_name, fl_name].iter().collect();
        let mut my_file = fs::File::create(file_path).expect("Error creating file");

    // pendulum parameters
    let period = 2.0;
    let omega = 2.0 * PI / period;
    let omega2 = omega * omega;

    // time parameters
    let tf = 100.0 * period;
    let dt_millis: u64 = 5;
    let dt = 0.001 * (dt_millis as f64);
    let mut delay_millis: u64;
    
    // initial conditions
    let mut t_0 = 0.0;
    let mut t_draw = 0.0;
    let mut t_start = 0.0;
    let mut t_end = 0.0;
    
    let mut theta_0 = 175.0 *  PI / 180.0;
    let mut phi_0 = 90.0 *  PI / 180.0;
    let theta_dot_0 = 0.0;
    let phi_dot_0 = 0.0;
    let mut p_0 = 8.0/3.0 * theta_dot_0 + phi_dot_0 * (theta_0 - phi_0).cos();
    let mut q_0 = 2.0/3.0 * phi_dot_0 + theta_dot_0 * (theta_0 - phi_0).cos();

    //total energy
    let mut h = right_hand_side(theta_0, phi_0, p_0, q_0, omega2)[4];

    // draw the pendulum initial state
    
    while !is_key_down(KeyCode::Space) {
        clear_background(LIGHTGRAY);

        let par = get_window_par(screen_height(), screen_width());

        draw_pendulum(theta_0, phi_0, par);

        let w = par[1];

        draw_text(format!("Use left/right or up/down arrows to move the pendulum!").as_str(), 5.0 * w, 5.0 * w, 3.0 * w, BLACK);
        draw_text(format!("Press SPACE to start!").as_str(), 5.0 * w, 10.0 * w, 3.0 * w, BLACK);
        
        let ft = get_frame_time() as f64;

        if !is_key_down(KeyCode::Left) {
            theta_0 = theta_0 + SPEED * ft;
        }

        if !is_key_down(KeyCode::Right) {
            theta_0 = theta_0 - SPEED * ft;
        }

        if !is_key_down(KeyCode::Down) {
            phi_0 = phi_0 + SPEED * ft;
        }

        if !is_key_down(KeyCode::Up) {
            phi_0 = phi_0 - SPEED * ft;
        }


        next_frame().await;
    }

    // introduce the state variables
    let mut theta: f64;
    let mut phi: f64;
    let mut p: f64;
    let mut q: f64;

    //start the calculation

    while t_0 < tf {

        delay_millis = (1000.0 * (t_end - t_start)) as u64;
        
        if dt_millis > delay_millis {thread::sleep(time::Duration::from_millis(dt_millis - delay_millis))};

        t_start = get_time();

        t_0 = t_0 + dt;

        if t_0 - t_draw > 0.016666667 {
            t_draw = t_0;

            clear_background(LIGHTGRAY);

            let par = get_window_par(screen_height(), screen_width());
            let w = par[1];

            draw_pendulum(theta_0, phi_0, par);

            draw_text(format!("Time elapsed, seconds: {:.3}", t_0).as_str(), 5.0 * w, 5.0 * w, 3.0 * w, BLACK);
            draw_text(format!("Total energy, arb. units: {:.5}", h).as_str(), 5.0 * w, 10.0 * w, 3.0 * w, BLACK);

            writeln!(my_file, "{} {} {} {}", t_0, angle_in_degrees(theta_0), angle_in_degrees(phi_0), h).expect("Error writing to file");
       
            next_frame().await;
        }
        

        let rk_1 = right_hand_side(theta_0, phi_0, p_0, q_0, omega2);
        theta = theta_0 + rk_1[0] * dt / 2.0;
        phi = phi_0 + rk_1[1] * dt / 2.0;
        p = p_0 + rk_1[2] * dt / 2.0;
        q = q_0 + rk_1[3] * dt / 2.0;

        let rk_2 = right_hand_side(theta, phi, p, q, omega2);
        theta = theta_0 + rk_2[0] * dt / 2.0;
        phi = phi_0 + rk_2[1] * dt / 2.0;
        p = p_0 + rk_2[2] * dt / 2.0;
        q = q_0 + rk_2[3] * dt / 2.0;

        let rk_3 = right_hand_side(theta, phi, p, q, omega2);
        theta = theta_0 + rk_3[0] * dt;
        phi = phi_0 + rk_3[1] * dt;
        p = p_0 + rk_3[2] * dt;
        q = q_0 + rk_3[3] * dt;

        let rk_4 = right_hand_side(theta, phi, p, q, omega2);

        theta_0 = theta_0 + (rk_1[0] + 2.0 * rk_2[0] + 2.0 * rk_3[0] + rk_4[0]) * dt / 6.0;
        phi_0 = phi_0 + (rk_1[1] + 2.0 * rk_2[1] + 2.0 * rk_3[1] + rk_4[1]) * dt / 6.0;
        p_0 = p_0 + (rk_1[2] + 2.0 * rk_2[2] + 2.0 * rk_3[2] + rk_4[2]) * dt / 6.0;
        q_0 = q_0 + (rk_1[3] + 2.0 * rk_2[3] + 2.0 * rk_3[3] + rk_4[3]) * dt / 6.0;

        h = rk_1[4];

        t_end = get_time();
            
        
    }
}