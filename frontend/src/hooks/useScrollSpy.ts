import { useEffect } from 'react';

interface UseScrollSpyOptions {
  containerRef: React.RefObject<HTMLElement>;
  selector: string;
  onIntersect: (id: string) => void;
  dependencies?: any[];
  rootMargin?: string;
  threshold?: number;
}

export function useScrollSpy({
  containerRef,
  selector,
  onIntersect,
  dependencies = [],
  rootMargin = '-40% 0px -40% 0px',
  threshold = 0,
}: UseScrollSpyOptions) {
  useEffect(() => {
    if (!containerRef.current) return;
    const observer = new IntersectionObserver(
      (entries) => {
        entries.forEach((entry) => {
          if (entry.isIntersecting) {
            const id = entry.target.getAttribute('data-id');
            if (id) onIntersect(id);
          }
        });
      },
      {
        root: containerRef.current,
        rootMargin,
        threshold,
      }
    );
    const elements = containerRef.current.querySelectorAll(selector);
    elements.forEach((el) => observer.observe(el));
    
    return () => observer.disconnect();
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [containerRef, selector, rootMargin, threshold, ...dependencies]);
}
