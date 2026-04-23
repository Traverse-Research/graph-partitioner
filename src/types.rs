use std::fmt;
use std::result::Result as StdResult;

pub type Idx = i32;
pub type Real = f32;
pub const NOPTIONS: usize = 40;
pub type Result<T> = StdResult<T, Error>;

#[derive(Debug, PartialEq, Eq)]
pub enum Error {
    Input,
    Memory,
    Other,
}

impl fmt::Display for Error {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Error::Input => write!(f, "input error"),
            Error::Memory => write!(f, "memory error"),
            Error::Other => write!(f, "unknown error"),
        }
    }
}

impl std::error::Error for Error {}

impl From<NewGraphError> for Error {
    fn from(_: NewGraphError) -> Self {
        Error::Input
    }
}

impl From<NewMeshError> for Error {
    fn from(_: NewMeshError) -> Self {
        Error::Input
    }
}

#[derive(Debug)]
pub struct InvalidGraphError {
    msg: &'static str,
}

impl fmt::Display for InvalidGraphError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.msg)
    }
}

#[derive(Debug)]
#[non_exhaustive]
pub enum NewGraphError {
    NoConstraints,
    NoParts,
    TooLarge,
    InvalidGraph(InvalidGraphError),
}

impl NewGraphError {
    pub(crate) fn msg(msg: &'static str) -> Self {
        NewGraphError::InvalidGraph(InvalidGraphError { msg })
    }
}

impl fmt::Display for NewGraphError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            NewGraphError::NoConstraints => write!(f, "ncon must be positive"),
            NewGraphError::NoParts => write!(f, "nparts must be positive"),
            NewGraphError::TooLarge => write!(f, "array length does not fit in Idx"),
            NewGraphError::InvalidGraph(e) => write!(f, "{}", e),
        }
    }
}

impl std::error::Error for NewGraphError {}

#[derive(Debug)]
pub struct InvalidMeshError {
    msg: &'static str,
}

impl fmt::Display for InvalidMeshError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.msg)
    }
}

#[derive(Debug)]
#[non_exhaustive]
pub enum NewMeshError {
    NoParts,
    TooLarge,
    InvalidMesh(InvalidMeshError),
}

impl NewMeshError {
    pub(crate) fn msg(msg: &'static str) -> Self {
        NewMeshError::InvalidMesh(InvalidMeshError { msg })
    }
}

impl fmt::Display for NewMeshError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            NewMeshError::NoParts => write!(f, "nparts must be positive"),
            NewMeshError::TooLarge => write!(f, "array length does not fit in Idx"),
            NewMeshError::InvalidMesh(e) => write!(f, "{}", e),
        }
    }
}

impl std::error::Error for NewMeshError {}
