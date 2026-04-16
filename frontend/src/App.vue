<template>
  <div class="app">
    <nav class="navbar">
      <div class="nav-brand">
        <span class="brand-icon">🧬</span>
        <span class="brand-name">MetagenApp</span>
      </div>
      <div class="nav-links">
        <router-link to="/">Nueva corrida</router-link>
        <router-link to="/history">Historial</router-link>
      </div>
    </nav>
    <main class="main-content">
      <router-view />
    </main>
  </div>
</template>

<style>
*, *::before, *::after { box-sizing: border-box; margin: 0; padding: 0; }

body {
  font-family: 'Inter', 'Segoe UI', system-ui, sans-serif;
  background: #0f1117;
  color: #e2e8f0;
  min-height: 100vh;
}

.app { display: flex; flex-direction: column; min-height: 100vh; }

.navbar {
  display: flex;
  align-items: center;
  justify-content: space-between;
  padding: 0 2rem;
  height: 60px;
  background: #1a1d27;
  border-bottom: 1px solid #2d3148;
  position: sticky;
  top: 0;
  z-index: 100;
}

.nav-brand {
  display: flex;
  align-items: center;
  gap: 0.5rem;
  font-size: 1.1rem;
  font-weight: 700;
  color: #7c83fd;
}

.brand-icon { font-size: 1.4rem; }

.nav-links { display: flex; gap: 1.5rem; }

.nav-links a {
  color: #94a3b8;
  text-decoration: none;
  font-size: 0.9rem;
  font-weight: 500;
  padding: 0.3rem 0.6rem;
  border-radius: 6px;
  transition: all 0.15s;
}

.nav-links a:hover,
.nav-links a.router-link-active {
  color: #7c83fd;
  background: rgba(124, 131, 253, 0.1);
}

.main-content {
  flex: 1;
  padding: 2rem;
  max-width: 960px;
  width: 100%;
  margin: 0 auto;
}
</style>
