import { createRouter, createWebHistory } from 'vue-router'
import NewRun    from '../views/NewRun.vue'
import Monitor   from '../views/Monitor.vue'
import Results   from '../views/Results.vue'
import History   from '../views/History.vue'

const routes = [
  { path: '/',                   component: NewRun },
  { path: '/monitor/:jobId',     component: Monitor, props: true },
  { path: '/results/:jobId',     component: Results, props: true },
  { path: '/history',            component: History },
]

export default createRouter({
  history: createWebHistory('/metagenapp'),
  routes,
})
