/* A small browser contract fixture for exercising the real static shell. */
'use strict';
const fs = require('node:fs');
const path = require('node:path');
const vm = require('node:vm');

function boot(options = {}) {
  const ids = new Map(), downloads = [], blobs = new Map(), frames = new Map(), timers = [];
  let nextFrame = 0, nextURL = 0;
  const decode = s => s.replace(/&amp;/g, '&').replace(/&lt;/g, '<').replace(/&gt;/g, '>').replace(/&quot;/g, '"');
  class Element {
    constructor(tag) { this.tagName = tag.toUpperCase(); this.attributes = new Map(); this.children = []; this.parentNode = null; this.listeners = new Map(); this._text = ''; this._value = null; this.checked = false; this.disabled = false; this.offsetWidth = 240; this.offsetHeight = 180; }
    setAttribute(name, value) {
      this.attributes.set(name, String(value));
      if (name === 'id') ids.set(String(value), this);
      if (name === 'checked') this.checked = true;
      if (name === 'disabled') this.disabled = true;
      if (name === 'value' && this.tagName === 'INPUT' && this._value === null) this._value = String(value);
    }
    getAttribute(name) { return this.attributes.has(name) ? this.attributes.get(name) : null; }
    removeAttribute(name) { this.attributes.delete(name); if (name === 'disabled') this.disabled = false; }
    get firstChild() { return this.children[0] || null; }
    get textContent() { return this._text + this.children.map(c => c.textContent).join(''); }
    set textContent(value) { this._text = String(value); this.children.forEach(c => { c.parentNode = null; }); this.children = []; }
    get value() {
      if (this.tagName === 'SELECT') {
        const opts = this.getElementsByTagName('option');
        if (this._value !== null) return opts.some(o => o.value === this._value) ? this._value : '';
        const option = opts.find(o => o.getAttribute('selected') !== null) || opts[0]; return option ? option.value : '';
      }
      if (this.tagName === 'OPTION') return this.getAttribute('value') || this.textContent;
      return this._value === null ? this.getAttribute('value') || '' : this._value;
    }
    set value(value) { this._value = String(value); }
    appendChild(child) { child.parentNode = this; this.children.push(child); return child; }
    removeChild(child) { const index = this.children.indexOf(child); if (index < 0) throw new Error('missing child'); this.children.splice(index, 1); child.parentNode = null; return child; }
    getElementsByTagName(tag) { const found = []; const target = tag.toUpperCase(); for (const child of this.children) { if (child.tagName === target || target === '*') found.push(child); found.push(...child.getElementsByTagName(tag)); } return found; }
    getElementsByClassName(name) { return this.getElementsByTagName('*').filter(e => (e.getAttribute('class') || '').split(/\s+/).includes(name)); }
    getBoundingClientRect() { return { left: 10, top: 20, width: 800, height: 420, right: 810, bottom: 440 }; }
    addEventListener(name, fn) { if (!this.listeners.has(name)) this.listeners.set(name, []); this.listeners.get(name).push(fn); }
    dispatch(name, detail = {}) { const event = { target: this, currentTarget: this, defaultPrevented: false, preventDefault() { this.defaultPrevented = true; }, ...detail }; for (const fn of this.listeners.get(name) || []) fn(event); return event; }
    focus() { if (document.activeElement && document.activeElement !== this) document.activeElement.dispatch('blur'); document.activeElement = this; this.dispatch('focus'); }
    click() { if (this.disabled) return; if (this.tagName === 'A') downloads.push({ filename: this.getAttribute('download'), href: this.getAttribute('href'), blob: blobs.get(this.getAttribute('href')) }); this.dispatch('click'); }
  }
  const document = { hidden: false, activeElement: null, listeners: new Map(), getElementById: id => ids.get(id) || null,
    createElement: tag => new Element(tag), createElementNS: (_, tag) => new Element(tag),
    getElementsByClassName: name => document.documentElement.getElementsByClassName(name),
    addEventListener(name, fn) { if (!this.listeners.has(name)) this.listeners.set(name, []); this.listeners.get(name).push(fn); },
    dispatch(name) { for (const fn of this.listeners.get(name) || []) fn({ target: document }); }
  };
  const shell = new Element('document'), stack = [shell], html = fs.readFileSync(path.join(__dirname, 'index.html'), 'utf8');
  const voidTags = new Set(['META', 'LINK', 'INPUT', 'BR', 'HR', 'IMG']);
  for (const token of html.matchAll(/<\/?([a-z][a-z0-9-]*)([^>]*?)>|([^<]+)/gi)) {
    if (token[3] !== undefined) { stack[stack.length - 1]._text += decode(token[3]); continue; }
    if (token[0].startsWith('</')) { if (stack.length > 1) stack.pop(); continue; }
    const el = new Element(token[1]);
    for (const attr of token[2].matchAll(/([^\s=/'">]+)(?:\s*=\s*(?:"([^"]*)"|'([^']*)'|([^\s>]+)))?/g)) el.setAttribute(attr[1], decode(attr[2] ?? attr[3] ?? attr[4] ?? ''));
    stack[stack.length - 1].appendChild(el); if (!voidTags.has(el.tagName) && !token[0].endsWith('/>')) stack.push(el);
  }
  document.documentElement = shell.getElementsByTagName('html')[0]; document.body = shell.getElementsByTagName('body')[0];
  const sandbox = { document, console, Blob, URL: { createObjectURL(blob) { const url = 'blob:fixture-' + ++nextURL; blobs.set(url, blob); return url; }, revokeObjectURL(url) { blobs.delete(url); } },
    matchMedia: () => ({ matches: !!options.reducedMotion }), setTimeout: fn => { timers.push(fn); return timers.length; },
    requestAnimationFrame: fn => { const id = ++nextFrame; frames.set(id, fn); return id; }, cancelAnimationFrame: id => frames.delete(id) };
  sandbox.window = sandbox; sandbox.globalThis = sandbox;
  vm.createContext(sandbox);
  for (const filename of ['field.js', 'fem.js', 'app.js']) vm.runInContext(fs.readFileSync(path.join(__dirname, filename), 'utf8'), sandbox, { filename });
  return { document, sandbox, ids, downloads, blobs, frames, timers, $: id => ids.get(id),
    input(id, value, event = 'input') { const el = ids.get(id); el.value = String(value); el.dispatch(event); },
    frame(time) { const callbacks = [...frames.values()]; frames.clear(); callbacks.forEach(fn => fn(time)); },
    runTimers() { timers.splice(0).forEach(fn => fn()); },
    elements() { return document.documentElement.getElementsByTagName('*'); } };
}
module.exports = { boot };
