---
layout: page
title: Physics
permalink: /physics/
---

## Pendulums

<ul>
{% for post in site.categories.pendulums %}
  <li>
    <a href="{{ post.url }}">{{ post.title }}</a>
  </li>
{% endfor %}
</ul>

<!-- ## Kinematics

<ul>
{% for post in site.categories.kinematics %}
  <li>
    <a href="{{ post.url }}">{{ post.title }}</a>
  </li>
{% endfor %}
</ul> -->
