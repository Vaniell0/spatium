// Primary module unit for spatium.core. Re-exports all partitions so consumers
// can `import spatium.core;` and receive the full surface.
export module spatium.core;
export import :concepts;
export import :error;
export import :epsilon;
export import :precision;
export import :verify;
