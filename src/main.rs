use fft::{fft, generate_wn};
use std::f64::consts::PI;  

pub const SAMPLE_LEN: usize = 16;
pub const SAMPLE_LEN_F: f64 = SAMPLE_LEN as f64;                  
pub const W: f64 = 2.0 * PI / SAMPLE_LEN_F;

// generate sinusoid with harmonic of sample period for creating t
// fundamental frequency is harmonic = 1                          
fn generate_sinusoid(harmonic: i64, mag: f64, phase: f64) -> Vec<f64> { 
    let harmonic_f: f64 = harmonic as f64;                        
    let mut output: Vec<f64> = vec![0.0; SAMPLE_LEN];        
    let mut n: f64;                                               
    let mut phi: f64;                                             
                                                                  
    for index in 0..SAMPLE_LEN {                                  
        n = index as f64;                                         
        phi = W * harmonic_f * n + phase;                         
        output[index] = mag * phi.cos();                          
    }                                                             
                                                                  
    output                                                
}  

// generate a test signal for FFT input                           
fn test_signal() -> Vec<f64> {                           
    let sig1: Vec<f64> = generate_sinusoid(1, 1.0, 0.0); 
    let sig2: Vec<f64> = generate_sinusoid(2, 0.5, PI);  
    let mut output: Vec<f64> = vec![0.0; SAMPLE_LEN];        

    for n in 0..SAMPLE_LEN {                                      
        output[n] = sig1[n] + sig2[n];                            
    }                                                             

    output                                                
}

fn main() {
    let wn: Vec<(f64, f64)> = generate_wn(SAMPLE_LEN);
    let input: Vec<f64> = test_signal();
    let output: Vec<(f64, f64)> = fft(input, wn, SAMPLE_LEN);
    println!("main output: {:?}", output);
}

